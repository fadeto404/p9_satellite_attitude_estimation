classdef MM_MEKF
    properties
        dx_pre {mustBeNumeric}      % Pre-update dx = [a', db']'
        dx_post {mustBeNumeric}     % Post-update dx = [a', db']'
        x_pre {mustBeNumeric}       % Pre-update full state vector
        x {mustBeNumeric}           % Post-reset full state vector x(+)
        omega_est {mustBeNumeric}   % Estimated ang. vel. for propagation
        q_ref {mustBeNumeric}       % Reference quaternion, postupdate q_ref is the estimate
        beta_g {mustBeNumeric}      % Gyro bias estimate
        a_g {mustBeNumeric}         % Gibbs vector representation of att. error
        dbeta_g {mustBeNumeric}     % Gyro bias error estimate
        A_att {mustBeNumeric}       % Attitude matrix

        Pxx_pre {mustBeNumeric}     % Pre-update covariance
        Pxx_post {mustBeNumeric}    % Post-update covariance

        H {mustBeNumeric}           % Approximated measurement transformation
        F {mustBeNumeric}           % Approximated process dynamics
        Q {mustBeNumeric}           % Process noise covariance matrix
        R {mustBeNumeric}           % Measurement noise covariance matrix
        G {mustBeNumeric}           % Some matrix ???
        K {mustBeNumeric}           % Kalman Gain Matrix

        z {mustBeNumeric}           % Most recent measurements
        omega {mustBeNumeric}       % Observed angular velocity wrt. body frame
        dt {mustBeNumeric}          % Time step magnitude [s]
        a_obs {mustBeNumeric}       % Observed attitude, Gibbs vector rep.
    end
    methods
        function obj=MM_MEKF(x0, P0, q0, omega0, F, H, G, Q, R, t_s)
            obj.dx_pre = x0; % x0 = [a_g0; b0]
            obj.Pxx_pre = P0; 
            obj.q_ref = q0;
            obj.F = F; % State propagation matrix, 6x6
            obj.H = H; % Measurement sensitivity matrix, 3x6
            obj.G = G; % Something something covariance?!
            obj.Q = Q; % Process covariance/spectral density? 6x6
            obj.R = R; % Measurement noise covariance/spectral density? 3x3
            obj.a_g = zeros(3,1);
            obj.dx_post = zeros(6,1); % [a_g; db_g] (Gibbs vector, gyro bias error)
            obj.Pxx_post = P0;
            obj.z = zeros(7,1); % [q_obs; omega]
            obj.omega = omega0;
            obj.dt = t_s; % Time step, [s]
            obj.beta_g = zeros(3,1); % TODO
            obj.dbeta_g = zeros(3,1);
        end
        function obj=propagate(obj)
            % Propagate reference quaternion using last gyro reading as
            % angular velocity
            obj.omega_est = obj.omega - obj.beta_g;
            norm_omega_est = norm(obj.omega_est);
            %q_ref_dot = sh_quat_mult(0.5*[obj.omega_est; 0], obj.q_ref);
            %obj.q_ref = obj.q_ref + q_ref_dot*obj.dt; % NO: see 6.60 for discrete quaternion propagation

            % Discretized quaternion propagation
            c = cos(0.5*norm_omega_est*obj.dt);
            psi = (sin(0.5*norm_omega_est*obj.dt)/norm_omega_est)*obj.omega_est;
            psi_cross = cross_prod_mat(psi);
            Theta = [c*eye(3)-psi_cross, psi; -psi', c];
            obj.q_ref = Theta*obj.q_ref;
            
            % Propagate x according to linearized dynamics
%             Phi = eye(6) + obj.F*obj.delta_t;
            Omega = cross_prod_mat(obj.omega_est);
            Omega_sqr = Omega^2;
            Phi11 = eye(3,3) - Omega*(sin(norm_omega_est*obj.dt)/norm_omega_est) + Omega_sqr * ((1-cos(norm_omega_est*obj.dt))/norm_omega_est^2);
            Phi12 = Omega * ((1-cos(norm_omega_est*obj.dt))/norm_omega_est^2) - eye(3)*obj.dt - Omega_sqr * ((norm_omega_est*obj.dt - sin(norm_omega_est*obj.dt))/norm_omega_est^3);
            Phi21 = zeros(3,3);
            Phi22 = eye(3,3);
            Phi = [Phi11, Phi12; Phi21, Phi22];
            % Try cond(Phi)/cond(F)
            %obj.dx_pre = Phi*obj.dx_post; % ???
            
            % Propagate covariance matrix according to Riccati equation
%             Pxx_dot = obj.F*obj.Pxx_post + obj.Pxx_post*obj.F' + obj.G*obj.Q*obj.G';
%             obj.Pxx_pre = obj.Pxx_post + Pxx_dot*obj.delta_t; % Discretize

            % Discrete-time covariance propagation
            sigma_v = sqrt(10)*10^(-7); % TODO: Gyro noise std. dev.
            sigma_u = sqrt(10)*10^(-10); % TODO: Gyro bias random walk std. dev.
            Q_11 = (sigma_v^2 * obj.dt + sigma_u^2 * obj.dt^3 / 3)*eye(3);
            Q_12 = -(sigma_u^2 * obj.dt^2 / 2)*eye(3);
            Q_22 = (sigma_u^2 * obj.dt)*eye(3);
            Q_k = [Q_11, Q_12; Q_12, Q_22]; % G*Q_k*G'
            obj.Pxx_pre = Phi*obj.Pxx_post*Phi' + obj.G*Q_k*obj.G';

            % Attitude matrix
            obj.A_att = quat_att_mat(obj.q_ref);

            % Compute Kalman gain, before new measurements arrive, saves
            % time (in a "real" implementation)
            H_a = obj.H(:,1:3);
            P_a = obj.Pxx_pre(1:3, 1:3);
            P_c = obj.Pxx_pre(1:3, 4:6);
            obj.K = [P_a; P_c] * H_a' / (H_a*P_a*H_a' + obj.R); 
        end
        function obj=update(obj, z)
            obj.z = z; % z = [q_obs; omega_obs]
            obj.omega = z(5:7,1);
            % Star tracker gives a quaternion estimate of attitude, assumed
            % in star craft frame, we convert to Gibbs vector rep. of error
            % quaternion
            obj.a_obs = quat2gvr(sh_quat_mult(obj.z(1:4), quat_inv(obj.q_ref)));
            
            P_a = obj.Pxx_pre(1:3, 1:3);
            P_c = obj.Pxx_pre(1:3, 4:6);
            
            obj.dx_post = obj.K*(obj.a_obs - obj.dx_pre(1:3));

            obj.Pxx_post = obj.Pxx_pre - obj.K*[P_a, P_c]; % TODO: Change to numerically stable Joseph form
        end

        function obj=reset(obj)
            % Reset gyro bias error estimate
            obj.beta_g = obj.beta_g + obj.dx_post(4:6);
            
            % Reset quaternion error estimate
            % Normalize after rotation, avoids accumulating error in norm,
            % only possible with Gibbs vector representation
            delta_q = [obj.a_g; 2];
            obj.q_ref = sh_quat_mult(delta_q, obj.q_ref);
            obj.q_ref = obj.q_ref./norm(obj.q_ref);
            obj.a_g = zeros(3,1);

            % Reset state vector & error state vector
            obj.dx_pre = zeros(6,1);
            obj.dx_post = zeros(6,1);
            obj.x = [obj.a_g; obj.beta_g];
        end
    end
end
classdef CAL_MEKF
    properties
        dx_pre {mustBeNumeric}      % Pre-update dx = [a', db', ds', dku', dkl']'
        dx_post {mustBeNumeric}     % Post-update dx = [a', db', ds', dku', dkl']'
        x {mustBeNumeric}           % Post-reset full state vector x(+)
        x_pre {mustBeNumeric}       % Pre-update full state vector x(-)
        omega_est {mustBeNumeric}   % Estimated ang. vel. for propagation
        q_ref {mustBeNumeric}       % Reference quaternion, postupdate q_ref is the estimate
        beta_g {mustBeNumeric}      % Gyro bias estimate
        a_obs {mustBeNumeric}       % Observed attitude error, Gibbs vector rep.
        s_g {mustBeNumeric}         % Gyro scaling factors
        k_U {mustBeNumeric}         % Gyro misalignment factor - upper
        k_L {mustBeNumeric}         % Gyro misalignment factor - lower
        A_att {mustBeNumeric}       % Attitude matrix
        S_est {mustBeNumeric}       % Scaling and misalignment factor estimate
        U_est {mustBeNumeric}       % Matrix for upper k
        L_est {mustBeNumeric}       % Matrix for lower k

        Pxx_pre {mustBeNumeric}     % Pre-update error state covariance
        Pxx_post {mustBeNumeric}    % Post-update error state covariance

        z {mustBeNumeric}           % Most recent measurements
        omega {mustBeNumeric}       % Observed angular velocity wrt. body frame

        H {mustBeNumeric}           % Approximated measurement transformation
        F {mustBeNumeric}           % Approximated process dynamics
        Q {mustBeNumeric}           % Process noise covariance matrix/spectral density
        R {mustBeNumeric}           % Measurement noise covariance matrix
        G {mustBeNumeric}           % Determines how noise enters into states
        K {mustBeNumeric}           % Kalman Gain Matrix
        Q_k {mustBeNumeric}         % Discretized spectral density
        Q_kk {mustBeNumeric}        % G*Q_k*G'
        Phi {mustBeNumeric}         % Discretized state propagation matrix

        dt {mustBeNumeric}          % Time step magnitude [s] (constant)
        n {mustBeNumeric}           % Number of states in dx

        cond_obs {mustBeNumeric}    % Condition number of observability matrix
    end
    methods
        function obj=CAL_MEKF(x0, P0, q0, omega0, F, H, G, Q, R, t_s)
            obj.dt = t_s; % Time step, [s]

            obj.dx_pre = zeros(15,1); % dx0 = [a_g0; b0] (Attitude error parametrized by Gibbs vector, gyro bias error)
            obj.dx_post = zeros(15,1);
            obj.Pxx_pre = P0;       % Delta x covariance matrix initial condition
            obj.Pxx_post = P0; 
            obj.q_ref = q0;         % Attitude estimate initial condition
            obj.beta_g = x0(4:6);  % Gyro bias estimate initial condition
            obj.s_g = x0(7:9);     % Initial scaling factor estimate
            obj.k_U = x0(10:12);   % Initial upper misalignment estimate
            obj.k_L = x0(13:15);   % Initial lower misalignment estimate

            obj.F = F; % State propagation matrix, 15x15
            obj.H = H; % Measurement sensitivity matrix, 3x15
            obj.G = G; % Something something covariance?!, 15x15
            obj.Q = Q; % Process covariance/spectral density, 15x15
            % obj.Q = blkdiag(var_v*eye(3), var_u*eye(3));
            obj.R = R; % Measurement noise covariance/spectral density? 3x3
            obj.Q_k = obj.dt*obj.G*obj.Q*obj.G'; % Discrete process noise covariance
            
            obj.z = zeros(7,1); % [q_obs; omega] Vector of measurements
            obj.omega = omega0; % Observed angular velocities
            obj.n = 15;
        end

        function obj=propagate(obj)
            % Propagate reference quaternion using last gyro reading as
            % angular velocity
            obj.omega_est = (eye(3) - obj.S_est)*(obj.omega - obj.beta_g);
            norm_omega_est = norm(obj.omega_est);
            
            % Discretize linearized dynamics of state dx
            Omega = cross_prod_mat(obj.omega_est);
            obj.F = [-Omega, -(eye(3) - obj.S_est), -diag(obj.omega - obj.beta_g), -obj.U_est, -obj.L_est; zeros(12,15)]; 
            obj.G = [blkdiag(-(eye(3) - obj.S_est), eye(3)); zeros(9,6)];
%             script_A = obj.dt*[-obj.F, obj.G*obj.Q*obj.G';
%                         zeros(15,15), obj.F'];
%             script_B = expm(script_A);
%             Phi = script_B(16:end, 16:end)';
%             obj.Q_k = Phi*script_B(1:15, 16:end);
            Phi = eye(15) + obj.F*obj.dt;
            obj.Q_k = obj.dt*obj.G*obj.Q*obj.G';
%             obj.cond_obs = cond(obsv(Phi, obj.H));

            % Discrete-time covariance propagation
            obj.Pxx_pre = Phi*obj.Pxx_post*Phi' + obj.Q_k;
            obj.Phi = Phi;

            % Discretized quaternion propagation
            c = cos(0.5*norm_omega_est*obj.dt);
            psi = (sin(0.5*norm_omega_est*obj.dt)/norm_omega_est)*obj.omega_est;
            psi_cross = cross_prod_mat(psi);
            Theta = [c*eye(3)-psi_cross, psi; -psi', c];
            obj.q_ref = Theta*obj.q_ref;
            
            % Attitude matrix - only necessary for vector measurements
            %obj.A_att = quat_att_mat(obj.q_ref);

            % "Propagate" x
            obj.x_pre = obj.x;
        end

        function obj=update(obj, z)
            % Compute Kalman gain
            obj.K = obj.Pxx_pre*obj.H' / (obj.H*obj.Pxx_pre*obj.H' + obj.R);

            obj.z = z; % z = [q_obs; omega_obs]
            obj.omega = z(5:7,1);

            % Star tracker gives a quaternion estimate of attitude, assumed
            % in spacecraft frame, we convert to Gibbs vector rep. of error
            % quaternion
            delta_q = sh_quat_mult(obj.z(1:4), quat_inv(obj.q_ref));
            obj.a_obs = quat2gvr(delta_q); % Attitude error between observation and prediction
            
            % Update error state
            obj.dx_post = obj.K*(obj.a_obs);
           
            % Update error state covariance (Joseph form for stability)
            obj.Pxx_post = (eye(15)-obj.K*obj.H)*obj.Pxx_pre*(eye(15)-obj.K*obj.H)'+obj.K*obj.R*obj.K';

        end

        function obj=reset(obj)
            % Reset gyro bias error estimate
            obj.beta_g = obj.beta_g + obj.dx_post(4:6);

            % Reset scaling and misalignment factors
            obj.s_g = obj.s_g + obj.dx_post(7:9);
            obj.k_U = obj.k_U + obj.dx_post(10:12);
            obj.k_L = obj.k_L + obj.dx_post(13:15);

                        % Update matrices
            obj.S_est = [obj.s_g(1) obj.k_U(1) obj.k_U(2);
                         obj.k_L(1) obj.s_g(2) obj.k_U(3);
                         obj.k_L(2) obj.k_L(3) obj.s_g(3)];
            obj.U_est = [obj.omega(2) - obj.beta_g(2), obj.omega(3) - obj.beta_g(3), 0;
                         0, 0, obj.omega(3) - obj.beta_g(3);
                         0, 0, 0];
            obj.L_est = [0, 0, 0;
                         obj.omega(1) - obj.beta_g(1), 0, 0;
                         0, obj.omega(1) - obj.beta_g(1), obj.omega(2) - obj.beta_g(2)];
            
            % Reset quaternion error estimate
            % Normalize after rotation, avoids accumulating error in norm,
            % only possible with Gibbs vector representation
            delta_q = [0.5*obj.dx_post(1:3); 1];
            obj.q_ref = sh_quat_mult(delta_q, obj.q_ref);
            obj.q_ref = obj.q_ref./norm(obj.q_ref);

            % Reset state vector to current est. & error state vector to 0
            obj.dx_pre = zeros(15,1);
            obj.dx_post = zeros(15,1);
            obj.x = [obj.dx_post(1:3); obj.beta_g; obj.s_g; obj.k_U; obj.k_L];
        end
    end
end
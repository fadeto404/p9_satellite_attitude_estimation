classdef MM_MEKF_vec
    properties
        dx_pre {mustBeNumeric}      % Pre-update dx = [a', db']'
        dx_post {mustBeNumeric}     % Post-update dx = [a', db']'
        x {mustBeNumeric}           % Post-reset full state vector x(+)
        omega_est {mustBeNumeric}   % Estimated ang. vel. for propagation
        q_ref {mustBeNumeric}       % Reference quaternion, postupdate q_ref is the estimate
        beta_g {mustBeNumeric}      % Gyro bias estimate
        a_obs {mustBeNumeric}       % Observed attitude error, Gibbs vector rep.
        A_att {mustBeNumeric}       % Attitude matrix

        Pxx_pre {mustBeNumeric}     % Pre-update error state covariance
        Pxx_post {mustBeNumeric}    % Post-update error state covariance

        nmax {mustBeNumeric}        %
        z {mustBeNumeric}           % Most recent measurements
        omega {mustBeNumeric}       % Observed angular velocity wrt. body frame
        b_meas {mustBeNumeric}      %
        b_meas_est {mustBeNumeric}  %

        H {mustBeNumeric}           % Approximated measurement transformation
        F {mustBeNumeric}           % Approximated process dynamics
        Q {mustBeNumeric}           % Process noise covariance matrix/spectral density
        R {mustBeNumeric}           % Measurement noise covariance matrix
        R_full {mustBeNumeric}        % All Measurement noise covariance matrix
        G {mustBeNumeric}           % Determines how noise enters into states
        K {mustBeNumeric}           % Kalman Gain Matrix
        Q_k {mustBeNumeric}         % Discretized spectral density
        Q_kk {mustBeNumeric}        % G*Q_k*G'

        dt {mustBeNumeric}          % Time step magnitude [s] (constant)

        cond_obs {mustBeNumeric}    % Condition number of observability matrix
    end
    methods
        function obj=MM_MEKF_vec(dx0, P0, q0, omega0, nmax, F, H, G, Q, R_full, t_s)
            obj.dx_pre = zeros(6,1); % dx0 = [a_g0; b0]
            obj.Pxx_pre = P0; 
            obj.q_ref = q0;
            obj.F = F; % State propagation matrix, 6x6
            obj.H = H; % Measurement sensitivity matrix, 3x6
            obj.G = G; % Something something covariance?!
            obj.Q = Q; % Process covariance/spectral density? 6x6
            obj.R = R_full; % Measurement noise covariance/spectral density? 3x3
            obj.R_full = R_full;
            obj.dx_post = zeros(6,1); % [a_g; db_g] (Gibbs vector, gyro bias error)
            obj.Pxx_post = P0;
            obj.nmax = nmax; % Maximum amount of stars observable at once
            obj.z = zeros(103,1); % [meas; omega]
            obj.omega = omega0;
            obj.dt = t_s; % Time step, [s]
            obj.beta_g = dx0(4:6,:);
            sigma_v = sqrt(10)*10^(-7); % TODO: Make this an input argument
            sigma_u = sqrt(10)*10^(-10); % TODO: Make this an input argument
            Q_11 = (sigma_v^2 * obj.dt + sigma_u^2 * obj.dt^3 / 3)*eye(3);
            Q_12 = (sigma_u^2 * obj.dt^2 / 2)*eye(3); % Negative?
            Q_22 = (sigma_u^2 * obj.dt)*eye(3);
            obj.Q_k = [Q_11, Q_12; Q_12, Q_22]; 
            obj.Q_kk = obj.G*obj.Q_k*obj.G'; % See literature
        end
        function obj=propagate(obj)
            % Propagate reference quaternion using last gyro reading as
            % angular velocity
            obj.omega_est = obj.omega - obj.beta_g;
            norm_omega_est = norm(obj.omega_est);
            %q_ref_dot = sh_quat_mult(0.5*[obj.omega_est; 0], obj.q_ref);
            %obj.q_ref = obj.q_ref + q_ref_dot*obj.dt; % NO: see 6.60 for discrete quaternion propagation
            
            % Discretized, linearized dynamics of state dx
            Omega = cross_prod_mat(obj.omega_est);
            Omega_sqr = Omega^2;
            Phi11 = eye(3,3) - Omega*(sin(norm_omega_est*obj.dt)/norm_omega_est) + Omega_sqr * ((1-cos(norm_omega_est*obj.dt))/norm_omega_est^2);
            Phi12 = Omega * ((1-cos(norm_omega_est*obj.dt))/norm_omega_est^2) - eye(3)*obj.dt - Omega_sqr * ((norm_omega_est*obj.dt - sin(norm_omega_est*obj.dt))/norm_omega_est^3);
            Phi21 = zeros(3,3);
            Phi22 = eye(3,3);
            Phi = [Phi11, Phi12; 
                   Phi21, Phi22];
%             obj.F = [-Omega, -diag(obj.omega - obj.beta_g); zeros(3,6)]; 
%             Phi = eye(6) + obj.F*obj.dt;
            if isempty(obj.H)
                obj.cond_obs = Inf;
            else
                obj.cond_obs = cond(obsv(Phi, obj.H));
            end
            
            % Propagate covariance matrix according to Riccati equation
%             Pxx_dot = obj.F*obj.Pxx_post + obj.Pxx_post*obj.F' + obj.G*obj.Q*obj.G';
%             obj.Pxx_pre = obj.Pxx_post + Pxx_dot*obj.delta_t; % Discretize

            % Discrete-time covariance propagation
            obj.Pxx_pre = Phi*obj.Pxx_post*Phi' + obj.Q_kk;

            % Propagate x according to linearized dynamics
            % Discretized quaternion propagation
            c = cos(0.5*norm_omega_est*obj.dt);
            psi = (sin(0.5*norm_omega_est*obj.dt)/norm_omega_est)*obj.omega_est;
            psi_cross = cross_prod_mat(psi);
            Theta = [c*eye(3)-psi_cross, psi; -psi', c];
            obj.q_ref = Theta*obj.q_ref;

            % Attitude matrix - only necessary for vector measurements
            %obj.A_att = quat_att_mat(obj.q_ref);

        end
        function obj=update(obj, z)
            obj.z = z; % z = [meas; omega_obs]

            n = 3*obj.nmax;
            b=obj.z(1:n,:); % True body measurements
            bm=obj.z(n+1:2*n,:); % Actual (noisy) body measurements
            im=obj.z(2*n+1:3*n,:); % True inertial measurements
            av=obj.z(3*n+1:3*n+obj.nmax,:); % Availible stars
            obj.omega = obj.z(3*n+obj.nmax+1:3*n+obj.nmax+3,:); % Gyroscope measurements

            obj.A_att = quat_att_mat(obj.q_ref);

            obj.R = [];
            obj.H = [];
            obj.b_meas = [];
            obj.b_meas_est = [];

            [row,col] = find(av==1);
            if isempty(row) == 1
                return;
            end

            for j = 1:length(row)
                bm_e(:,j)=obj.A_att*im(row(j)*3-2:row(j)*3,:); % Estimated body measurement
                bm_e_cpm = cross_prod_mat(bm_e(:,j)); % Estimated body measurement cross product matrix

                obj.b_meas(j*3-2:j*3,1)=bm(row(j)*3-2:row(j)*3,:); % Availible body measurements
                obj.b_meas_est(j*3-2:j*3,1)=bm_e(:,j); % Availible estimated body measurements

                obj.H(j*3-2:j*3,:)=[bm_e_cpm zeros(3)]; % Output sensitivity matrix
                obj.R(j*3-2:j*3,j*3-2:j*3)=obj.R_full(row(j)*3-2:row(j)*3,row(j)*3-2:row(j)*3); % Availible measurement covariance matrix
            end

            % Compute Kalman gain
            obj.K = obj.Pxx_pre*obj.H' / (obj.H*obj.Pxx_pre*obj.H' + obj.R);

        %             % Update error state
            obj.dx_post = obj.K*(obj.b_meas - obj.b_meas_est);

            % Star tracker gives a quaternion estimate of attitude, assumed
            % in spacecraft frame, we convert to Gibbs vector rep. of error
            % quaternion
%             delta_q = sh_quat_mult(obj.z(1:4), quat_inv(obj.q_ref));
%             obj.a_obs = quat2gvr(delta_q); % Attitude error between observation and prediction
%             
%             % Update error state
%             obj.dx_post = obj.K*(obj.a_obs - obj.dx_pre(1:3));
           
            % Update error state covariance (Joseph form for stability)
            obj.Pxx_post = (eye(6)-obj.K*obj.H)*obj.Pxx_pre*(eye(6)-obj.K*obj.H)'+obj.K*obj.R*obj.K';
        end

        function obj=reset(obj)
            % Reset gyro bias error estimate
            obj.beta_g = obj.beta_g + obj.dx_post(4:6);
            
            % Reset quaternion error estimate
            % Normalize after rotation, avoids accumulating error in norm,
            % only possible with Gibbs vector representation
            delta_q = [0.5*obj.dx_post(1:3); 1];
            obj.q_ref = sh_quat_mult(delta_q, obj.q_ref);
            obj.q_ref = obj.q_ref./norm(obj.q_ref);

            % Reset state vector to current est. & error state vector to 0
            obj.dx_pre = zeros(6,1);
            obj.dx_post = zeros(6,1);
            obj.x = [obj.q_ref; obj.beta_g];
        end
    end
end
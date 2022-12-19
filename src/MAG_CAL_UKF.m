classdef MAG_CAL_UKF < handle
    properties
        x_pre {mustBeNumeric}   % Pre-update state estimate
        x_post {mustBeNumeric}  % Post-update state estimate
        P_pre {mustBeNumeric}   % Pre-update error covariance
        P_post {mustBeNumeric}  % Post-update error covariance
        
        Q {mustBeNumeric}   % Process noise covariance
        Sig {mustBeNumeric} % Magnetometer measurement noise covariance
        R {mustBeNumeric}   % Scalar measurement noise covariance
        K {mustBeNumeric}   % Kalman gain matrix

%         f {function_handle} % (Non-)linear state propagation function
%         h {function_handle} % (Non-)linear measurement function

        alpha {mustBeNumeric}   % Point spread parameter
        kappa {mustBeNumeric}   % Point spread parameter 2
        beta {mustBeNumeric}    % Certainty-of-Gaussian parameter
    end % properties

    methods
        function obj = MAG_CAL_UKF(x0, P0, Sigma, alpha, kappa, beta)
            obj.x_pre = x0;
            obj.P_pre = P0;
            obj.Sig = Sigma;
            obj.alpha = alpha;
            obj.kappa = kappa;
            obj.beta = beta;
        end % constructor

        function [x, P] = propagate(obj)
            obj.x_pre = obj.x_post;
            obj.P_pre = obj.P_post;
            x = obj.x_pre;
            P = obj.P_pre;
        end % propagate

        function [x, P] = update(obj, Bk, Rk)
            obj.R = obj.compute_sigma_sqr(obj.x_pre, Bk, obj.Sig);
            z = norm(Bk)^2 - norm(Rk)^2; % Synthetic scalar measurement
            [z_hat, P_zz, P_xz] = obj.unscented_transform(obj.x_pre, obj.P_pre, @obj.h, obj.alpha, obj.kappa, obj.beta, Bk);
            obj.K = P_xz / (P_zz+obj.R);

            obj.x_post = obj.x_pre + obj.K*(z-z_hat);
%             obj.P_post = obj.P_pre - obj.K*P_xz';
            obj.P_post = obj.P_pre - obj.K*(P_zz+obj.R)*obj.K';

            x = obj.x_post;
            P = obj.P_post;
        end % update
    end % methods

    methods(Static)
        % Unscented transform: For estimating the mean of a stochastic variable
        % propagated by a nonlinear function g(x)
        function [mean_sig, sig_cov, Pxy] = unscented_transform(x, P, g, alpha, kappa, beta, B)
            % Setup
            n = size(x, 1); % No. of states
            sigma_points = zeros(n, 2 * n + 1); % Preallocate space for sigma points
            W = zeros(1, 2 * n + 1); % Weights
            lambda = alpha^2 * (n + kappa) - n; % Parameter for transform
%             [U,Lambda,~]    = svd(P);
%             C               = U*sqrt(Lambda);
            C = chol(P); % Compute Cholesky factor of covariance matrix
            C = chol(C*C');
        
            % Compute sigma points
            sigma_points(:, 1) = x; % Set first sigma point as the mean
            W(1) = lambda / (n + lambda); % Set first weight
            WP = W(1) + 1 - alpha^2 + beta;    
            for i = 1:n
                sigma_points(:, i + 1) = x + sqrt(n + lambda) * C(:, i); % Positive sigma point
                sigma_points(:, n + i + 1) = x - sqrt(n + lambda) * C(:, i); % Negative sigma point
                W(i + 1) = 1 / (2 * (n + lambda));
                W(n + i + 1) = W(i + 1);
            end
            
            % Transform sigma points
            f_sigma_points = zeros(1, 2*n+1);
            for i=1:2*n+1
                f_sigma_points(:,i) = g(sigma_points(:,i), B);
            end
            
            % Compute mean and covariance
            mean_sig = zeros(1,1);
            for i=1:2*n+1
                mean_sig = mean_sig + W(i)*f_sigma_points(:,i);
            end
            sig_cov = WP * (f_sigma_points(:, 1) - mean_sig) * (f_sigma_points(:, 1) - mean_sig)';
            for i = 2:2*n+1
                sig_cov = sig_cov + W(i) * (f_sigma_points(:, i) - mean_sig) * (f_sigma_points(:, i) - mean_sig)';
            end
            
            % Compute P_xy
            Pxy = WP * (sigma_points(:, 1) - x) * (f_sigma_points(:, 1) - mean_sig)';
            for i=2:2*n+1
                Pxy = Pxy + W(i) * (sigma_points(:, i) - x) * (f_sigma_points(:, i) - mean_sig)';
            end
        end % unscented_transform

        function [z_hat] = h(x, B)
            % Measurement model for the scalar measurement, evaluated at
            % last state estimate x_k with magnetometer measurement B_{k+1}
            S = [B(1)^2, B(2)^2, B(3)^2, 2*B(1)*B(2), 2*B(1)*B(3), 2*B(2)*B(3)]';
            D_mat = [x(4), x(7), x(8);
                     x(7), x(5), x(9);
                     x(8), x(9), x(6)];
            E_mat = 2*D_mat + D_mat^2;
            E_vec = [E_mat(1,1), E_mat(2,2), E_mat(3,3), E_mat(1,2), E_mat(1,3), E_mat(2,3)]';
            z_hat = [2*B', -S']*[(eye(3)+D_mat)*x(1:3); E_vec] - norm(x(1:3))^2;
        end % h (measurement function)

        function [var] = compute_sigma_sqr(x, B, Sigma)
            D_mat = [x(4), x(7), x(8);
                     x(7), x(5), x(9);
                     x(8), x(9), x(6)];
            var = 4*(((eye(3)+D_mat)*B - x(1:3))'*Sigma*((eye(3)+D_mat)*B - x(1:3))) + 2*trace(Sigma^2);
        end % compute_sigma_sqr

    end % methods(Static)
end % classdef

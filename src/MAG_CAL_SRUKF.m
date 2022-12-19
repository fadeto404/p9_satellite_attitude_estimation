classdef MAG_CAL_SRUKF < handle
    properties
        x_pre {mustBeNumeric}   % Pre-update state estimate
        x_post {mustBeNumeric}  % Post-update state estimate
        P_post {mustBeNumeric}  % Post-update error covariance
        S_pre {mustBeNumeric}   % Pre-update covariance Cholesky factor
        S_post {mustBeNumeric}  % Post-update covariance Cholesky factor

        X {mustBeNumeric}       % Propagated sigma points
        
        Q {mustBeNumeric}   % Process noise covariance
        Sig {mustBeNumeric} % Magnetometer measurement noise covariance
        R {mustBeNumeric}   % Scalar measurement noise covariance
        K {mustBeNumeric}   % Kalman gain matrix
        sqrR {mustBeNumeric}    % (Matrix) square root of R

        alpha {mustBeNumeric}   % Point spread parameter
        kappa {mustBeNumeric}   % Point spread parameter 2
        beta {mustBeNumeric}    % Certainty-of-Gaussian parameter
        lambda {mustBeNumeric}  % Intermediate parameter
        c {mustBeNumeric}       % Constant parameter for sigma point generation
        Wy {mustBeNumeric}      % Weights for mean in UT
        Wc {mustBeNumeric}      % Weights for covariance in UT
    end % properties

    methods
        function obj = MAG_CAL_SRUKF(x0, P0, Sigma, alpha, kappa, beta)
            obj.x_post = x0;
            obj.P_post = P0;
            obj.S_post = chol(obj.P_post);
            obj.Sig = Sigma; % Not sigma points, Sigma=magnetometer noise covariance
            obj.alpha = alpha;  % Typically: 0.1
            obj.kappa = kappa;  % 3-n (n=no. of states)
            obj.beta = beta;    % Typically 2 --> Assumes a Gaussian dist.
            n = length(x0);     % No. of states
            obj.lambda = obj.alpha^2*(n+obj.kappa)-n;
            obj.Wy = zeros(1, 2 * n + 1);
            obj.Wy(1) = obj.lambda / (n + obj.lambda);
            obj.Wy(2:end) = 1 / (2 * (n + obj.lambda));
            obj.Wc = obj.Wy;
            obj.Wc(1) = obj.Wc(1) + 1 - obj.alpha^2 + obj.beta;
            obj.c = sqrt(n+obj.lambda);
            obj.Q = zeros(9,9); % Process noise covariance, constant x -> 0
        end % constructor

        function [x, S] = update(obj, Bk, Rk)
            z = Bk(1)^2 + Bk(2)^2 + Bk(3)^2 - (Rk(1)^2 + Rk(2)^2 + Rk(3)^2); % Synthetic scalar measurement, Bk'*Bk - Rk'*Rk
            obj.R = obj.compute_sigma_sqr(obj.x_pre, Bk, obj.Sig); % Variance of synthetic measurement
            obj.sqrR = sqrtm(obj.R); % Std. dev. of synthetic measurement
            [z_hat, S_z, Z] = obj.unscented_transform(obj.x_pre, obj.X, @obj.h, 1, obj.sqrR, obj.Wy, obj.Wc, Bk);
            
            % Compute P_xz
            n = length(obj.x_pre);
%             P_xz = (obj.X-obj.x_pre)*diag(obj.Wc)*(Z-z_hat)';
            P_xz = obj.Wc(1) * (obj.X(:, 1) - obj.x_pre) * (Z(:, 1) - z_hat)';
            for i=2:2*n+1
                P_xz = P_xz + obj.Wc(i) * (obj.X(:, i) - obj.x_pre) * (Z(:, i) - z_hat)';
            end
            
            obj.K = (P_xz/S_z)/S_z';

            obj.x_post = obj.x_pre + obj.K*(z-z_hat);
            U = obj.K*S_z;
            obj.S_post = cholupdate(obj.S_pre, U, '-');

            x = obj.x_post;
            S = obj.S_post;
        end % update

        function [x_hat, S_x] = propagate(obj)
            chi = obj.generate_sigma_points(obj.x_post, obj.S_post, obj.c);
            n = length(obj.x_post);
            [x_hat, S_x, Xi] = obj.unscented_transform(obj.x_post, chi, @obj.f, n, obj.Q, obj.Wy, obj.Wc, zeros(3,1));
            obj.X = Xi; 
            obj.x_pre = x_hat;
            obj.S_pre = S_x;
        end
    end % methods

    methods(Static)

        function [sigma_points] = generate_sigma_points(x, S, c)
            n = length(x);

            A = c*S'; % Matrix with scaled rows of S
            Y = x(:,ones(1,n)); % Array of [x, x, x, ..., x] of size nxn
            sigma_points = [x Y+A Y-A]; 
            
%             sigma_points = zeros(n, 2 * n + 1); % Preallocate space for sigma points
%             sigma_points(:, 1) = x; % Set first sigma point as the mean
%             for i = 1:n
%                 sigma_points(:, i + 1) = x + c*S(:, i)'; % Positive sigma point
%                 sigma_points(:, n + i + 1) = x - c*S(:, i)'; % Negative sigma point
%             end

        end

        % Unscented transform: For estimating the mean of a stochastic variable
        % propagated by a nonlinear function g(x)
        function [mean_sig, sig_cov, Y] = unscented_transform(x, sigma_points, g, m, sqrR, Wy, Wc, B)
            % Setup
            n = length(x); % No. of states
            m = m; % No. of function outputs; g: R^n -> R^m
            
            % Transform sigma points
            f_sigma_points = zeros(m, 2*n+1);
            for i=1:2*n+1
                f_sigma_points(:,i) = g(sigma_points(:,i), B);
            end
            
            % Compute mean and covariance
            mean_sig = zeros(m,1);
            for i=1:2*n+1
                mean_sig = mean_sig + Wy(i)*f_sigma_points(:,i);
            end
            Y = f_sigma_points;
            
            S_y = qr([sqrt(Wc(2))*(f_sigma_points(:,2:end) - mean_sig), sqrR]', 0);
            sig_cov = cholupdate(S_y, sqrt(abs(Wc(1)))*(f_sigma_points(:,1) - mean_sig), '+');            
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
%             z_hat = [2*B', -S']*[(eye(3)+D_mat)*x(1:3); E_vec] - norm(x(1:3))^2;
            z_hat = -S'*E_vec + 2*B'*(eye(3) + D_mat)*x(1:3) - (x(1)^2 + x(2)^2 + x(3)^2);
        end % h (measurement function)

        function  [x_hat] = f(x, ~)
            x_hat = x;
        end

        function [var] = compute_sigma_sqr(x, B, Sigma)
            D_mat = [x(4), x(7), x(8);
                     x(7), x(5), x(9);
                     x(8), x(9), x(6)];
            var = 4*(((eye(3)+D_mat)*B - x(1:3))'*Sigma*((eye(3)+D_mat)*B - x(1:3))) + 2*trace(Sigma^2);
        end % compute_sigma_sqr
    end % methods(Static)
end % classdef
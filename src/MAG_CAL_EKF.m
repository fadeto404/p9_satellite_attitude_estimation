classdef MAG_CAL_EKF < handle
    properties
        x_post      % x = [b1, b2, b3, D11, D22, D33, D12, D13, D23]'
        P_post      % Error covariance matrix
        
        z           % Synthetic scalar measurement (not actual magnetometer output)
        z_hat       % Expected scalar measurement
        sigma_sqr   % Synthetic scalar measurement covariance
        
        H           % Measurement sensitivity matrix
        R           % Noise covariance matrix of magnetometer measurements
        K           % Kalman gain matrix
    end

    methods
        function obj = MAG_CAL_EKF(x0, P0, R)
            obj.x_post = x0;
            obj.P_post = P0;
            obj.R = R;
        end

        function [x, P] = update(obj, Bk, Rk)
            obj.z = Bk(1)^2 + Bk(2)^2 + Bk(3)^2 - Rk(1)^2 - Rk(2)^2 - Rk(3)^2; % Synthetic scalar measurement
            obj.H = obj.compute_H(obj.x_post, Bk); % Sensitivity matrix using B_k and x_{k-1}
            obj.z_hat = obj.h(obj.x_post, Bk); % Estimated scalar measurement using B_k and x_{k-1}
            obj.sigma_sqr = obj.compute_sigma_sqr(obj.x_post, Bk, obj.R); % Compute current measurement covariance

            obj.K = obj.P_post*obj.H'/(obj.H*obj.P_post*obj.H' + obj.sigma_sqr); % Kalman gain
%             obj.P_post = (eye(9) - obj.K*obj.H)*obj.P_post; % Update error covariance
            obj.P_post = (eye(9)-obj.K*obj.H)*obj.P_post*(eye(9)-obj.K*obj.H)' + obj.K*obj.sigma_sqr*obj.K';
            obj.x_post = obj.x_post + obj.K*(obj.z - obj.z_hat); % Update state estimate
            x = obj.x_post;
            P = obj.P_post;
        end

        function [x, P] = propagate(obj)
            % Given that x_dot = 0, propagation does not change anything
            % (the estimated quantities are constants)
            x = obj.x_post;
            P = obj.P_post;
        end
    end
    methods(Static)
        function [z_hat] = h(x, B)
            % Measurement model for the scalar measurement, evaluated at
            % last state estimate x_k with magnetometer measurement B_{k+1}
            S = [B(1)^2, B(2)^2, B(3)^2, 2*B(1)*B(2), 2*B(1)*B(3), 2*B(2)*B(3)]';
            D_mat = [x(4), x(7), x(8);
                     x(7), x(5), x(9);
                     x(8), x(9), x(6)];
            D_vec = x(4:9);
            E_mat = 2*D_mat + D_mat^2;
            E_vec = [E_mat(1,1), E_mat(2,2), E_mat(3,3), E_mat(1,2), E_mat(1,3), E_mat(2,3)]';
            z_hat = -S'*E_vec + 2*B'*(eye(3) + D_mat)*x(1:3) - (x(1)^2 + x(2)^2 + x(3)^2);
%             z_hat = -B'*E_mat*B + 2*B'*(eye(3) + D_mat)*x(1:3) - norm(x(1:3))^2;
%             z_hat = [2*B', -S']*[(eye(3)+D_mat)*x(1:3); E_vec] - norm(x(1:3))^2;
        end
        function [H] = compute_H(x, B)
            % Compute sensitivity matrix given state estimate x_k and
            % magnetometer measurement B_{k+1}
            S = [B(1)^2, B(2)^2, B(3)^2, 2*B(1)*B(2), 2*B(1)*B(3), 2*B(2)*B(3)]';
            D_mat = [x(4), x(7), x(8);
                     x(7), x(5), x(9);
                     x(8), x(9), x(6)];
            M11 = diag(2*x(4:6));
            M21 = [x(7), x(7), 0;
                   x(8), 0, x(8);
                   0, x(9), x(9)];
            M12 = 2*M21';
            M22 = [x(4)+x(5), x(9), x(8);
                   x(9), x(4)+x(6), x(7);
                   x(8), x(7), x(5)+x(6)];
            M_ED = 2*eye(6) + [M11, M12; M21, M22];
            J = [B(1)*x(1), B(2)*x(2), B(3)*x(3), B(1)*x(2)+B(2)*x(1), B(1)*x(3)+B(3)*x(1), B(2)*x(3)+B(3)*x(2)];
            H = [2*B'*(eye(3) + D_mat)-2*x(1:3)', -S'*M_ED + 2*J];
        end
        function [var] = compute_sigma_sqr(x, B, Sigma)
            D_mat = [x(4), x(7), x(8);
                     x(7), x(5), x(9);
                     x(8), x(9), x(6)];
            var = 4*(((eye(3)+D_mat)*B - x(1:3))'*Sigma*((eye(3)+D_mat)*B - x(1:3))) + 2*trace(Sigma^2);
        end
    end
end
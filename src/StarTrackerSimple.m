classdef StarTrackerSimple
    properties
        axis_noise_var {mustBeNumeric} % Variance 3x1; [x,y,z] in ST frame
        body_frame_tf {mustBeNumeric} % Rotation mat. from ST to BF
        cov_mat {mustBeNumeric} % Noise cov. mat. in BF
        cov_mat_chol_factor {mustBeNumeric} % Cholesky factor of cov. mat.
    end
    methods
        function obj = StarTrackerSimple(noise_var_arr, bf_rot)
            obj.axis_noise_var = noise_var_arr;
            obj.body_frame_tf = bf_rot;
            obj.cov_mat = diag(noise_var_arr);
            obj.cov_mat_chol_factor = chol(obj.cov_mat)';
        end
        function [q_meas, q_err] = simulate_reading(self,q_true)            
            % Noise is modelled as Gaussian dist. random angles with
            % Generate noise angles (gamma -> Rx, beta -> Ry, alpha -> Rz)
            %ang_err = randn([1,3]).*sqrt(self.axis_noise_var);
            
            % Rotate errors into body frame
            %ang_err = self.body_frame_tf*ang_err;

            % Generate errors in star tracker frame
            ang_err = self.cov_mat_chol_factor*randn([3,1]);
            ang_err = self.body_frame_tf*ang_err;
            
            % Extrinsic rotation matrix for measurement error
            gamma = ang_err(1); % Error angle around body frame x
            beta = ang_err(2); % Error angle around body frame y
            alpha = ang_err(3); % Error angle around body frame z
            ca = cos(alpha); % Shorthands for sine and cosine
            sa = sin(alpha);
            cb = cos(beta);
            sb = sin(beta);
            cy = cos(gamma);
            sy = sin(gamma);
            R_err = [ca*cb, ca*sb*sy-sa*cy, ca*sb*cy+sa*sy;
                     sa*cb, sa*sb*sy+ca*cy, sa*sb*cy-ca*sy;
                     -sb, cb*sy, cb*cy];

            % Rotate into body frame
            %R_err = self.body_frame_tf*R_err;

            % Convert measurement noise rotation matrix to quaternion
            q_err = [R_err(2,3) - R_err(3,2);
                     R_err(3,1) - R_err(1,3);
                     R_err(1,2) - R_err(2,1);
                     1 + trace(R_err)];
            q_err = q_err./norm(q_err);
            
            % Rotate true quaternion by the measurement noise quaternion
            q_meas = sh_quat_mult(q_err, q_true);
            q_meas = q_meas./norm(q_meas);
        end
    end
end
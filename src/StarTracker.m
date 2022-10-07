classdef StarTracker
    properties
        axis_noise {mustBeNumeric}
        %body_frame_transform {mustBeNumeric}
    end
    methods
        function obj = StarTracker(noise_arr)
            obj.axis_noise = noise_arr;
            %body_frame_transform = transform_mat;
        end
        function [q_meas, q_err] = simulate_reading(self,q_true)
            %q_true = quaternion(randn(),randn(),randn(),randn()).normalize();
            q_true_eul = quat2eul(q_true, 'ZYX');
            %variances = [deg2rad(5/(60^2)), deg2rad(5/(60^2)), deg2rad(50/(60^2))];
            q_eul_noise = randn([1,3]).*sqrt(self.axis_noise);
            q_meas_eul = q_true_eul + q_eul_noise;
            q_meas = quaternion(q_meas_eul,'euler', 'ZYX', 'frame').normalize();
            q_err = q_true - q_meas;
        end
    end
end
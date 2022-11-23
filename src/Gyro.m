classdef Gyro
    properties
        bias {mustBeNumeric}
        bias_covariance {mustBeNumeric}
        noise_covariance {mustBeNumeric}
        delta_t {mustBeNumeric}
        bias_cov_chol {mustBeNumeric}
        noise_cov_chol {mustBeNumeric}
    end
    methods
        function obj = Gyro(noise_cov, bias_cov, b0, dt)
            obj.bias = b0;
            obj.bias_covariance = bias_cov;
            obj.noise_covariance = noise_cov;
            obj.delta_t = dt;
            obj.bias_cov_chol = chol(obj.bias_covariance)';
            obj.noise_cov_chol = chol(obj.noise_covariance)';
        end
        function [omega_meas, bias, obj] = simulate_reading(obj, omega_true)
            obj = obj.propagate_bias();
            % w_out = w_true + eta_v + eta_u, 
            % eta_v ~ N(0,R_v), eta_u(t) ~ N(µ_u,R_u), µ_u = b0
            omega_meas = omega_true + obj.noise_cov_chol*randn(3,1)*obj.delta_t + obj.bias*obj.delta_t;
            bias = obj.bias;
        end
        function obj = propagate_bias(obj)
            obj.bias = obj.bias + (obj.bias_cov_chol*randn(3,1))*obj.delta_t;
        end
    end
end
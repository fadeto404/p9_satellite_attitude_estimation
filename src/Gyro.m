classdef Gyro
    % See Fundamentals of Spacecraft Attitude Determination and 
    % Control, by Markley & Crassidis 2014, pp. 147
    properties
        bias {mustBeNumeric} % True bias
        bias_var {mustBeNumeric} % Continuous time bias variance
        noise_var {mustBeNumeric} % Continuous time noise variance
        dt {mustBeNumeric} % Discrete time step
        disc_noise {mustBeNumeric} % Discretized noise variance
        disc_bias {mustBeNumeric} % Discretized bias variance
        S {mustBeNumeric} % Scale factor and misalignment matrix
        IS {mustBeNumeric} % I_3 + S
    end
    methods
        function obj = Gyro(noise_var, bias_var, b0, dt, S)
            if nargin < 5
                S = zeros(3,3); % No scaling and/or misalignment
            end
            obj.bias = b0;
            obj.dt = dt;
            obj.S = S;
            obj.IS = eye(3) + obj.S;
            obj.bias_var = bias_var(1);
            obj.noise_var = noise_var(1);
            obj.disc_noise = sqrt(obj.noise_var/obj.dt + obj.bias_var*obj.dt/12);
            obj.disc_bias = sqrt(obj.bias_var*obj.dt);
        end

        function [omega_meas, bias, obj] = simulate_reading(obj, omega_true)
            bias_prev = obj.bias;
            obj = obj.propagate_bias();
            bias = obj.bias;
            % w_out = (I+S)*w_true + eta_v + eta_u, 
            % eta_v ~ N(0,R_v), eta_u(t) ~ N(µ_u,R_u), µ_u = b0
            omega_meas = obj.IS*omega_true + ...
                         obj.disc_noise*randn(3,1) + ...
                         0.5*(bias_prev + bias);
        end

        function obj = propagate_bias(obj)
            obj.bias = obj.bias + obj.disc_bias*randn(3,1);
        end
    end
end
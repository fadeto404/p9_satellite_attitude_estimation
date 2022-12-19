classdef Magnetometer < handle
    properties
        D_true {mustBeNumeric}      % Symmetric matrix of scaling factors and non-orthogonality corrections
        b_true {mustBeNumeric}      % True bias
        OT {mustBeNumeric}          % Rotation matrix containing info on axis misalignments
        noise_cov {mustBeNumeric}   % Noise covariance matrix, discrete-time
        noise_cov_chol {mustBeNumeric}
        IpDinv {mustBeNumeric}      % (I_3 + D_true)^(-1); For improved computation time
        B_k {mustBeNumeric}         % Most recent measurement
        eps_k {mustBeNumeric}       % Most recent noise vector
        dt  {mustBeNumeric}         % Discretization time step
    end
    methods
        function obj = Magnetometer(D_true, b_true, OT, Sigma_k)
            obj.D_true = D_true;
            obj.b_true = b_true;
            obj.OT = OT;
            obj.noise_cov = Sigma_k;
            obj.noise_cov_chol = chol(Sigma_k)';
            obj.IpDinv = inv(eye(3) + D_true);
        end

        function [B_k] = simulate_reading(obj, A_true, R_k)
            obj.eps_k = obj.noise_cov_chol*randn(3,1);
            obj.B_k = obj.IpDinv*(obj.OT*A_true*R_k + obj.b_true + obj.eps_k);
            B_k = obj.B_k;
        end
    end
end
classdef QuaternionFuser
    % QuaternionFuser Uses weighted quaternion average to fuse n quaternions
    %   Fuse a set of quaternions of arbitrary size, using weighted
    %   quaternion averaging. For unweighted average set R as n identity
    %   matrices, e.g. R(:,:,i) = eye(3) for i=1,...,n, n=4. Using covariance
    %   matrices for R results in Maximum Likelihood Estimate of the 
    %   average quaternion. Quaternions are expected as column vectors of 4
    %   doubles, arranged in an array such that q(1:4, i) has corresponding
    %   weighting matrix R(1:3,1:3,i)
    properties
        R_bar {mustBeNumeric} % Approximate fused quaternion covariance
        R_inv {mustBeNumeric} % Array of inverted weighting matrices
        q_bar {mustBeNumeric} % Estimated quaternion
    end
    methods
        function obj=QuaternionFuser(R)
            obj.set_weighting_matrices(R);
            obj.q_bar = zeros(4,1);
        end

        function obj=fuse(obj, q_array)
            M = zeros(4,4);
            [~, num_quats] = size(q_array);
            for i=1:num_quats
                Xi = quat_xi_mat(q_array(:,i));
                M = M - Xi*obj.R_inv(:,:,i)*Xi';
            end
            [eigvec, eigval] = eig(M);
            [~, ind] = max(diag(eigval));
            obj.q_bar = eigvec(:,ind);
        end

        function obj=set_weighting_matrices(obj, R)
            [n, m, num_quats] = size(R);
            obj.R_inv = zeros(n,m,num_quats);
            for i=1:num_quats
                obj.R_inv(:,:,i) = inv(R(:,:,i));
            end
            obj.R_bar = inv(sum(obj.R_inv, 3));
        end
    end
end
classdef QuaternionFuser
    % QuaternionFuser Uses weighted quaternion average to fuse n quaternions
    %   Fuse a set of quaternions of arbitrary size, using weighted
    %   quaternion averaging. For unweighted average set R as n identity
    %   matrices, e.g. R(:,:,1:4) = eye(3) for n=4. Using covariance
    %   matrices for R results in Maximum Likelihood Estimate of the 
    %   average quaternion.
    properties
        R_bar {mustBeNumeric}
        R_inv {mustBeNumeric}
        q_bar {mustBeNumeric}
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
            [~, ~, num_quats] = size(R);
            obj.R_inv = zeros(3,3,num_quats);
            for i=1:num_quats
                obj.R_inv(:,:,i) = inv(R(:,:,i));
            end
            obj.R_bar = inv(sum(obj.R_inv, 3));
        end
    end
end
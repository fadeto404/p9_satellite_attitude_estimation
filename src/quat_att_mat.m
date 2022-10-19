function A = quat_att_mat(q)
% Convert quaternion to attitude matrix
% Input: 4x1 array with quaternion coefficients, q(1:3): the vector part,
% q(4): the scalar part.
% Mathematical formulation:
% A = (q(4)^2 - norm(q)^2)*eye(3) - 2*q(4)*cross_prod_mat(q(1:3)) + 2*q(1:3)*q(1:3)';

    A = [q(1)^2-q(2)^2-q(3)^2+q(4)^2, 2*(q(1)*q(2) + q(3)*q(4)), 2*(q(1)*q(3) - q(2)*q(4));
         2*(q(1)*q(2) - q(3)*q(4)), -q(1)^2+q(2)^2-q(3)^2+q(4)^2, 2*(q(3)*q(2) + q(1)*q(4));
         2*(q(1)*q(3) + q(2)*q(4)), 2*(q(3)*q(2) - q(1)*q(4)), -q(1)^2-q(2)^2+q(3)^2+q(4)^2];
end
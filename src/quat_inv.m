% Quaternion inversion, q must be 4x1, q(1:3) vector part, q(4) scalar part
function q_inv = quat_inv(q)
    q_inv = [-q(1:3); q(4)]./norm(q);
end
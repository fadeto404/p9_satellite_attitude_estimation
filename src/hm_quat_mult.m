% Quaternion multiplication by Hamilton's convention
function res = hm_quat_mult(q, p)
    res = [p(4)*q(1:3) + q(4)*p(1:3) + cross(q(1:3), p(1:3));
           q(4)*p(4) - dot(q(1:3), p(1:3))];
end
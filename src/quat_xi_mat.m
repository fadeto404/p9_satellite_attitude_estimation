% Matrix part of Hamilton quaternion product matrix (vector part = q)
function Xi = quat_xi_mat(q)
    Xi = [q(4), -q(3), q(2);
             q(3), q(4), -q(1);
             -q(2), q(1), q(4);
             -q(1), -q(2), -q(3)];
end
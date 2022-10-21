% Matrix part of Shuster quaternion product matrix (vector part = q)
function Psi = quat_psi_mat(q)
    Psi = [q(4), q(3), -q(2);
           -q(3), q(4), q(1);
           q(2), -q(1), q(4);
           -q(1), -q(2), -q(3)];
end
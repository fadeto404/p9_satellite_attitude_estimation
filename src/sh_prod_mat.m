% Quaternion product matrix by Shuster's multiplication convention
function sh_mat = sh_prod_mat(q)
    sh_mat = [quat_psi_mat(q), q];
end
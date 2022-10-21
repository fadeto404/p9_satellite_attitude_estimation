% Quaternion product matrix by Hamilton's multiplication convention
function hm_mat = hm_prod_mat(q)
    hm_mat = [quat_xi_mat(q), q];
end
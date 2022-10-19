% Compute the cross product matrix of 3-component column vector x
function xc = cross_prod_mat(x)
    xc = zeros([3,3]);
    xc(1,2) = -x(3);
    xc(2,1) = x(3);
    xc(3,1) = -x(2);
    xc(1,3) = x(2);
    xc(3,2) = x(1);
    xc(2,3) = -x(1);
end
% Convert quaternion to Gibbs vector representation (a_g = 2*g = 2*q/q4)
function a_g=quat2gvr(q)
    a_g = 2*q(1:3)./q(4);
end
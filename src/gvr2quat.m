% Convert Gibbs vector representation to quaternion, normalization optional
function q=gvr2quat(a_g, normalize=0)
    if normalize==0
        q = [a_g; 2];
    else
        q = [a_g; 2] * (1/sqrt(4 + a_g(1)^2 + a_g(2)^2 + a_g(3)^2));
    end
end
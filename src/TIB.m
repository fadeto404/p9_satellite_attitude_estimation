function out = TIB(phi,theta,psi)
ct = cos(theta);
st = sin(theta);
cp = cos(phi);
sp = sin(phi);
cs = cos(psi);
ss = sin(psi);

out = [ct*cs, sp*st*cs-cp*ss, cp*st*cs+sp*ss;
       ct*ss, sp*st*ss+cp*cs, cp*st*ss-sp*cs;
       -st, sp*ct,cp*ct];
end


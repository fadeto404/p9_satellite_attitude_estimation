function R = ext_rot_mat_xyz(gamma, beta, alpha)
    % Shorthands for sine and cosine
    cy = cos(gamma);
    sy = sin(gamma);
    cb = cos(beta);
    sb = sin(beta);
    ca = cos(alpha); 
    sa = sin(alpha);
    
    % Intermediate rotation matrices
    Rx = [1, 0, 0;
          0, cy, -sy;
          0, sy, cy];
    Ry = [cb, 0, sb;
          0 , 1, 0;
          -sb, 0, cb];
    Rz = [ca, -sa, 0;
          sa, ca, 0;
          0, 0, 1];
    R = Rz*Ry*Rx;
    
    % Resulting rotation matrix
%     R = [ca*cb, ca*sb*sy-sa*cy, ca*sb*cy+sa*sy;
%          sa*cb, -sa*sb*sy+ca*cy, -sa*sb*cy-ca*sy;
%          -sb, cb*sy, cb*cy];
end
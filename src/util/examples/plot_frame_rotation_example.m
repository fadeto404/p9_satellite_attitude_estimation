clear; clc; close all;

% Angles for extrinsic rotation, order: X-Y-Z
gamma = deg2rad(-90); % Angle around x
beta = deg2rad(0);   % Angle around y
alpha = deg2rad(0);  % Angle around z

% Shorthands for sine and cosine
cy = cos(gamma);
sy = sin(gamma);
cb = cos(beta);
sb = sin(beta);
ca = cos(alpha); 
sa = sin(alpha);

% Intermediate rotation matrices
% Rx = [1, 0, 0;
%       0, cy, -sy;
%       0, sy, cy];
% Ry = [cb, 0, sb;
%       0 , 1, 0;
%       -sb, 0, cb];
% Rz = [ca, -sa, 0;
%       sa, ca, 0;
%       0, 0, 1];
% R1 = Rz*Ry*Rx;

% Resulting rotation matrix
R1 = [ca*cb, ca*sb*sy-sa*cy, ca*sb*cy+sa*sy;
         sa*cb, sa*sb*sy+ca*cy, sa*sb*cy-ca*sy;
         -sb, cb*sy, cb*cy];

R1 = [1, 0, 0; 0, 0, 1; 0, 1, 0];

% Optional rotation matrix describing base frame relative to world frame
R0 = [1 0 0; 0 0 -1; 0 1 0];

plot_frame_rotation(R1) % Rotate base frame = world frame by R1
%plot_frame_rotation(R1, R0) % Rotate base frame = R0*world frame by R1
%plot_frame_rotation(R1*R0) % Rotate world frame by R1*R0

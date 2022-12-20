clear; clc; close all;

% Angles for extrinsic rotation, order: X-Y-Z
gamma = deg2rad(24); % Angle around x
beta = deg2rad(64);   % Angle around y
alpha = deg2rad(13);  % Angle around z

R1 = ext_rot_mat_xyz(gamma, beta, alpha);

% R1 = [1, 0, 0; 0, 0, 1; 0, 1, 0];

% Optional rotation matrix describing base frame relative to world frame
R0 = [1 0 0; 0 0 -1; 0 1 0];

plot_frame_rotation(R1) % Rotate base frame = world frame by R1
%plot_frame_rotation(R1, R0) % Rotate base frame = R0*world frame by R1
%plot_frame_rotation(R1*R0) % Rotate world frame by R1*R0

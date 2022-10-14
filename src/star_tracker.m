%% Star tracker simulation
clear; clc; close all;
% Assume stars in catalog are stationary, equidistant to spacecraft/earth, 
% uniformly distributed on a sphere

%% Generate star catalog: Uniform random points on a sphere
% Muller method
rng(2);
num_stars = 7542;
u = randn([num_stars, 1]);
v = randn([num_stars, 1]);
w = randn([num_stars, 1]);

unit = symunit;
r = 25*unit.ly; % Distance to stars: 25 lightyears
[r,~] = separateUnits(unitConvert(r,unit.km));
r = double(r);

P_norm = sqrt(u.*u + v.*v + w.*w);
P = r*[u./P_norm, v./P_norm, w./P_norm];

figure(1); hold on;
scatter3(P(:,1), P(:,2), P(:,3), Marker=".")
grid on; grid minor;
axis equal;

%% Given true star tracker attitude, pick stars (simulate physical system)
% Accuracy: a few arcseconds along boresight, more for rots. about
% boresight
% Update rate: 0.5-10Hz
% Geometry: Pinhole camera
% opt_ax = randn([3,1]);
% opt_ax = r*(opt_ax./norm(opt_ax))
% norm(opt_ax)
% plot3([0; opt_ax(1)], [0; opt_ax(2)], [0; opt_ax(3)], "LineWidth", 2)

% Camera is defocused --> stars cover several pixels (2x2 or 3x3) -->
% centroid is determined with accuracy kappa_cent
% kappa_cent = 0.1; % [px]
% f = 50E-6; % [mm] to [km]; Distance to focal plane
% res = 1024; % [px]; Resolution/image width and height, assuming square img
% u0 = res/2;
% v0 = res/2;
% img_coord_func = @(s)([u0 + f*s(1)/s(3); v0 + f*s(2)/s(3)]);


%% Estimate star tracker attitude (simulate measurements)
% Approach 1: Given true state, add slight noise to attitude measurements
x_true = quaternion(randn(),randn(),randn(),randn()).normalize();
% Rotation order: Z-Y-X
variances = [deg2rad(60/(60^2)), deg2rad(5/(60^2)), deg2rad(5/(60^2))];

st = StarTracker(variances)
%[x_meas, meas_err] = st.simulate_reading(x_true)
x_meas_hist = [quaternion(0,0,0,0)];
meas_err_hist = [quaternion(0,0,0,0)];
for i=1:200
    [x_meas_hist(i,:), meas_err_hist(i,:)] = st.simulate_reading(x_true);
    x_meas_eul(i,:) = quat2eul(x_meas_hist(i,:), 'ZYX')
    x_meas_eul(i,:) = x_meas_eul(i,:)./norm(x_meas_eul(i,:))
end

x_true_eul = quat2eul(x_true, 'ZYX');
%x_meas_eul = quat2eul(x_meas_hist, 'ZYX');
%x_noise = randn([1,3]).*sqrt(variances);
%x_meas_eul = x_true_eul + x_noise;
%x_meas = quaternion(x_meas_eul,'euler', 'ZYX', 'frame').normalize()
%meas_err = x_true - x_meas

opt_ax = r*(x_true_eul./norm(x_true_eul))
opt_ax_est = r*(x_meas_eul)
plot3([0; opt_ax(1)], [0; opt_ax(2)], [0; opt_ax(3)], "LineWidth", 2)
scatter3([opt_ax_est(:,1)], [opt_ax_est(:,2)], [opt_ax_est(:,3)], "LineWidth", 1, "Marker", '.')



























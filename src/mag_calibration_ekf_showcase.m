clc; clear; close all;

%% Orbit constants
ME = 5.97217e+24;   % Mass of earth [kg]
G =  6.67430e-11;   % Gravitational constant [m³/(kg.s²)]
GM = G*ME;          % Earth's: gravitational constant * mass [m³/s²]
RE = 6371.2*1000;   % Earth's radius: average used in IGRF13 [m]
nanoTesla = 400*1000;       % Satellite orbit radius = RE+r [m]
v0 = sqrt(GM/(RE+nanoTesla));   % Necessary speed to keep satellite in circular orbit at distance RE+r from center of earth
T_o = 2*pi*(RE+nanoTesla)/v0;   % Orbital period [s]
inclination = deg2rad(30); % Orbit inclination (angle of orbit relative to equator)
%% Simulation time setup
rng(42);
t_s = 0.1; % [s]
num_orbits = 1; % [.]
num_min = num_orbits*T_o/60; % Sim time in [min]
N = round(num_min*60 / t_s, 0); % No. of time steps [.]
T = 0:t_s:(N-1)*t_s; % [s]

%% Attitude simulation
% Generate angular accelerations [x,y,z]
omega_dot = 10*[0.000001*sin(2.3*pi*T./(60*60)); 0.000001*sin(2*pi*T./(60*60)); 0.000001*cos(2*pi*T./(60*60))];

% Generate angular velocities by Euler integration
omega = zeros(3,N);
for i=1:N-1
    omega(:,i+1) = omega(:,i) + t_s*omega_dot(:,i);
end

% CONSTANT VELOCITY
% w=0.1*pi/180*[sin(t_s*0.01*T') sin(t_s*0.0085*T') cos(t_s*0.0085*T')];
% omega = w';
% omega = ones(3,N).*[inv(91.5/(2*pi)*60)*0; -inv(91.5/(2*pi)*60); 0];

% True initial attitude quaternion
q0 = 1/sqrt(2)*[1; 0; 0; 1];
q_true = [q0, zeros(4,N-1)];
q_norm = [1, zeros(1,N-1)];

% Integrating quaternions
% Convert ang. vel. to unit quaternion, multiply prev. q_true by q_omega
for i=1:N-1
    %A = 0.5*sh_prod_mat([omega(:,i); 0]);
    if omega(:,i) == zeros(3,1)
        q_true(:,i+1) = q_true(:,i);
    else
        omega_norm = norm(omega(:,i));
        omega_normalized = omega(:,i)./omega_norm;
        q_omega = [omega_normalized*sin(0.5*omega_norm*t_s); cos(0.5*omega_norm*t_s)];
        q_true(:,i+1) = sh_quat_mult(q_omega, q_true(:,i));
    end
    % Compute norm and renormalize to compensate for effect of
    % computational inaccuracies
    q_norm(:,i+1) = norm(q_true(:,i+1));
    q_true(:,i+1) = q_true(:,i+1)/q_norm(:,i+1);
    q_norm(:,i+1) = norm(q_true(:,i+1));
end

%% Two-body problem simulation
pos0 = [RE+nanoTesla; 0; 0]; % Start position at intersection of equator and IERS reference meridian (WGS84) at distance RE+r from center of earth
vel0 = ext_rot_mat_xyz(inclination,0,0)*[0; v0; 0]; % Start velocity, a vector of magnitude v0 (necessary speed to keep circular orbit) perpendicular to gravitational force, rotated by inclination
orbit_state0 = [pos0; vel0];  % x, y, z, x_dot, y_dot, z_dot
[~,orbit_state] = ode89(@satellite_earth_twobp, T, orbit_state0);
orbit_state=orbit_state';

%% IGRM simulation
orbit_pos_geod = zeros(3,N); % (lat, lon, alt) [deg,deg,m]
R_mag = zeros(3,N); % (x_NED, y_NED, z_NED) [nT, nT, nT]
wgs84 = wgs84Ellipsoid("meter"); % Reference ellipsoid for geodetic coordinates
timestamp = datenum([2022 12 13 14 15 16]); % Timestamp for use in IGRF; TODO: propagate timestamp with simulation time
tic
for i=1:N
    orbit_pos_geod(:,i) = ecef2geodetic(wgs84, orbit_state(1,i),orbit_state(2,i),orbit_state(3,i))';    % Coordinate change
    R_mag(:,i) = igrfmagm(orbit_pos_geod(3,i), orbit_pos_geod(1,i), orbit_pos_geod(2,i), decyear(datetime('today')))'; % Simulate magnetic field according to IGRF-13
%     R_mag(:,i) = igrf(timestamp, orbit_pos_geod(1,i), orbit_pos_geod(2,i), orbit_pos_geod(3,i), 'geod');
end
toc

%% Convert measurements from NED to ECEF
% [x_e,y_e,z_e]=ned_to_ecef(R_mag(1,1), R_mag(2,1), R_mag(3,1), orbit_pos_geod(1,i), orbit_pos_geod(2,i), orbit_pos_geod(3,i))
% R_mag(:,1)
% for i=1:N
%     R_mag(:,i) = ned_to_ecef(R_mag(1,i), R_mag(2,i), R_mag(3,i), orbit_pos_geod(1,i), orbit_pos_geod(2,i), orbit_pos_geod(3,i));
% end

%% Magnetometer setup
unit = symunit;
milliGauss = 1*unit.mG;
[nanoTesla,~] = separateUnits(unitConvert(milliGauss,unit.nT));
nanoTesla = double(nanoTesla); % My dumb ass can't convert from mG to nT

D1 = [0.0800, 0.0520, 0.0500;
      0.0520, 0.0500, 0.0490;
      0.0500, 0.0490, 0.0750]; % Unitless, scaling factors!
b1 = [5.00; 6.00; 5.50]*nanoTesla;
x_true = [b1; D1(1,1); D1(2,2); D1(3,3); D1(1,2); D1(1,3); D1(2,3)];
OT1 = eye(3); % ext_rot_mat_xyz(deg2rad(0.05), deg2rad(0.02), deg2rad(0.01));
Sigma1 = eye(3).*((0.05*nanoTesla*ones(3,1)).^2);
mgm1 = Magnetometer(D1, b1, OT1, Sigma1);

%% Magnetometer simulation
B_mag = zeros(3,N);
for i=1:N
   B_mag(:,i) = mgm1.simulate_reading(quat_att_mat(q_true(:,i)), R_mag(:,i));
end

%% Calibration EKF setup
x0 = zeros(9,1);
% x0 = [b1; D1(1,1); D1(2,2); D1(3,3); D1(1,2); D1(1,3); D1(2,3)];
P0 = eye(9).*[1000; 1000; 1000; 0.01; 0.01; 0.01; 0.01; 0.01; 0.01]*1e-2;
% P0 = eye(9).*[10; 10; 10; 0.001; 0.01; 0.01; 0.01; 0.01; 0.01]*1e-2;
R_mgm1 = Sigma1;
mag_cal_ekf = MAG_CAL_EKF(x0, P0, R_mgm1);

%% Calibration UKF setup
% f = @(x)(0);
% h = mag_cal_ekf.h;
% Q_ukf = zeros(length(x0));
alpha = 0.1;
kappa = -6;
beta = 2;
mag_cal_ukf = MAG_CAL_UKF(x0, P0, R_mgm1, alpha, kappa, beta);

%% Calibration
x_sim = zeros(9,N);
P_sim = zeros(9,9,N);
for i=1:N
    [x_sim(:,i), P_sim(:,:,i)] = mag_cal_ekf.update(B_mag(:,i), R_mag(:,i));
%     [x_sim(:,i), P_sim(:,:,i)] = mag_cal_ukf.update(B_mag(:,i), R_mag(:,i));
%     [~,~] = mag_cal_ukf.propagate();
end

figure(1); 
subplot(3,1,1); hold on;
plot(T, x_sim(1:3,:))
plot(T, ones(size(x_sim(1:3,:))).*x_true(1:3), "LineStyle","--")
hold off;
subplot(3,1,2); hold on;
plot(T, x_sim(4:6,:))
plot(T, ones(size(x_sim(4:6,:))).*x_true(4:6), "LineStyle","--")
hold off;
subplot(3,1,3); hold on;
plot(T, x_sim(7:9,:))
plot(T, ones(size(x_sim(7:9,:))).*x_true(7:9), "LineStyle","--")
hold off;

b_hat = x_sim(1:3,end)
D_hat = [x_sim(4,end), x_sim(7,end), x_sim(8,end);
         x_sim(7,end), x_sim(5,end), x_sim(9,end);
         x_sim(8,end), x_sim(9,end), x_sim(6,end)]

err = x_sim(:,end) - x_true;
deviation = err./x_true * 100;
% mat2str(b_hat, 10)
% mat2str(deviation, 10)

%% Compare corrected measurements to true (not correct)
% D1_hat = [x_sim(4,end), x_sim(7,end), x_sim(8,end);
%          x_sim(7,end), x_sim(5,end), x_sim(9,end);
%          x_sim(8,end), x_sim(9,end), x_sim(6,end)];
% b1_hat = [x_sim(1,end); x_sim(2,end); x_sim(3,end)];
% B_mag_cal = zeros(3,N);
% for i=1:N
%     B_mag_cal(:,i) = (eye(3) + D1_hat)*B_mag(:,i) - b1_hat;
% end
% B_err = R_mag - B_mag_cal;
% meas_err = R_mag - B_mag;
% 
% figure(2); hold on;
% subplot(2,1,1);
% plot(T, B_err)
% subplot(2,1,2);
% plot(T, meas_err)
% hold off;
% 
% cal_mse = mean(B_err.^2, 2)
% raw_mse = mean(meas_err.^2, 2)


%% Animated plot
% Setting up the Plot
figure(3); hold on
title(sprintf('Trajectory\nTime: %0.0f sec', T(1)), 'Interpreter', 'Latex');
xlabel('x', 'Interpreter', 'Latex')
ylabel('y', 'Interpreter', 'Latex')
zlabel('z', 'Interpreter', 'Latex')
grid minor  % Adding grid lines
axis equal  % Equal axis aspect ratio
% view(-37.5,30);  % Setting viewing angle
load('topo.mat', 'topo', 'topolegend');
axesm('globe', 'Geoid', RE);
meshm(topo, topolegend); 
demcmap(topo);

% Create file name variable
% filename = 'animation.gif';
% Plotting with no color to set axis limits
plot3(orbit_state(1,:), orbit_state(2,:), orbit_state(3,:), 'Color','none');
% Plotting the first iteration
p = plot3(orbit_state(1,1),orbit_state(2,1),orbit_state(3,1),'b');
m = scatter3(orbit_state(1,1),orbit_state(2,1),orbit_state(3,1),'filled','k');
% Iterating through the length of the time array
for k = 1:300:N
    % Updating the line
    p.XData = orbit_state(1,1:k);
    p.YData = orbit_state(2,1:k);
    p.ZData = orbit_state(3,1:k);
    % Updating the point
    m.XData = orbit_state(1,k); 
    m.YData = orbit_state(2,k);
    m.ZData = orbit_state(3,k);
    % Updating the title
    title(sprintf('Trajectory\nTime: %0.2f sec', T(k)),...
    'Interpreter','Latex');
    % Delay
    pause(0.01)
    % Saving the figure
%     frame = getframe(gcf);
%     im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);
%     if k == 1
%         imwrite(imind,cm,filename,'gif', 'Loopcount',inf,...
%         'DelayTime',0.1);
%     else
%         imwrite(imind,cm,filename,'gif','WriteMode','append',...
%         'DelayTime',0.1);
%     end
end

%% Animation (from Sid, kills laptops)
% x_initial = orbit_state0(1);
% y_initial = orbit_state0(2);
% z_initial = orbit_state0(3);
% a = scatter3(x_initial, y_initial, z_initial,'MarkerEdgeColor','r','LineWidth',1,Marker="o");
% [lat,long,rad] = cart2sph(x_initial,z_initial,y_initial);
% [XYZ,H,D,I,F,DXDYDZ,DH,DD,DI,DF] = igrfmagm(rad,lat,long,decyear(datetime('today')));
% MagNorm = XYZ/norm(XYZ);
% b = plot3([0,MagNorm(1)*1e6],[0,MagNorm(2)*1e6],[0,MagNorm(3)*1e6],'color','k','linewidth',2);
% for i = T
%     x = nanoTesla*cos(i);
%     y = nanoTesla*sin(i)*cos(deg2rad(theta));
%     z = -nanoTesla*sin(i)*sin(deg2rad(theta)); 
%     drawnow;
%     rho = norm([x y z]);
%     thetaE = acos(z/rho);
%     phiE = 0;
%     psiE = atan2(y,x);
%     lat = 90-thetaE*180/pi;
%     long = psiE*180/pi;
%     alt = (rho - EarthRad);
% %     [lat,long,rad] = cart2sph(x,y,z);
%     [XYZ,H,D,I,F,DXDYDZ,DH,DD,DI,DF] = igrfmagm(alt,lat,long,decyear(datetime('today')));
%     BNED = [XYZ(1);XYZ(2);-XYZ(3)];
%     BI = TIB(phiE,thetaE+pi,psiE)*BNED;
%     MagNorm = BI/norm(BI);
%     delete(a);
%     delete(b);
%     a = scatter3(x, y, z,'MarkerEdgeColor','r','LineWidth',1,Marker="o");
%     b = plot3([x,x+MagNorm(1)*1e6],[y,y+MagNorm(2)*1e6],[z,z+MagNorm(3)*1e6],'color','k','linewidth',2);
% end
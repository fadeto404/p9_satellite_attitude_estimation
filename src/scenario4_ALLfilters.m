close all; clear; clc;

% All Sensors


%% MEKF showcase simulation with fused star tracker measurement and gyro
rng(42);
% Simulation time basics
t_s = 0.1; % [s]
num_min = 100; % Sim time in [min]
N = num_min*60 / t_s; % No. of time steps
T = 0:t_s:(N-1)*t_s;

%% Gyro setup
b0 = ones(3,1)*(0.1 / (180/pi*3600)); % 0.1 [deg/hr] to [rad/s] (init bias)
gyro_noise_variance = (sqrt(10)*10^(-7))^2; % Noise variance
gyro_bias_variance = (sqrt(10)*10^(-10))^2; % Bias variance
gyro_noise_cov = eye(3)*gyro_noise_variance; % Noise covariance matrix
gyro_bias_cov = eye(3)*gyro_bias_variance; % Bias random walk cov. mat.
gy1 = Gyro(gyro_noise_variance, gyro_bias_variance, b0, t_s); % Gyro sim init

%% Star tracker setup
% Vector Star Tracker Setup
sig_st=87.2665/3*1e-6; % Noise variance

%% Magnetometer sensor setup
sig_mag = 0.15^2; % nano Tesla

%% Fine Sun Sensor setup
sig_sun = deg2rad(2); % 0.5 - 2 degrees

%% Filter setup
nsensors = 3; % Number of sensors
nmax=nsensors*10; % Maximum amount of simultaneously observations

R_full=blkdiag(sig_mag^2*eye(3*(nmax/nsensors)),sig_sun^2*eye(3*(nmax/nsensors)),sig_st^2*eye(3*(nmax/nsensors))); % Full measurement covariance matrix 
% R_full=blkdiag(sig_mag^2*eye(3*(nmax/nsensors)),sig_sun^2*eye(3*(nmax/nsensors))); % Full measurement covariance matrix 
% R_full=blkdiag(sig_sun^2*eye(3*(nmax/nsensors))); % Full measurement covariance matrix 

dx0 = [zeros(3,1); zeros(3,1)];
P0 = [(0.1*pi/180)^2*eye(3), zeros(3); zeros(3), (0.2*pi/180/3600)^2*eye(3)];
P0 = blkdiag((0.1*pi/180/3600)^2*eye(3) * eye(3), (0.2*pi/180/3600)^2 * eye(3));

qest0 = 1/sqrt(2)*[1; 0; 0; 1];

d = 10^-10;
variances = [deg2rad(d), deg2rad(d), deg2rad(d)];
st1 = StarTrackerSimple(variances, eye(3));
q_est0 = st1.simulate_reading(qest0);

omega0 = 0.0000001*ones(3,1);
% omega0 = [inv(91.5/(2*pi)*60)*0; -inv(91.5/(2*pi)*60); 0]; % For constant velocity sim
% omega0 = 0.1*pi/180*[sin(t_s*0.01*T(1)') sin(t_s*0.0085*T(1)') cos(t_s*0.0085*T(1)')]';
F0 = [zeros(3) -eye(3); zeros(3,6)];
G = [-eye(3) zeros(3); zeros(3) eye(3)];
H = [eye(3) zeros(3)];
R = R_full;
Q = [gyro_noise_cov, zeros(3); zeros(3), gyro_bias_cov];

%% MEKF Object

mekf_obj = MM_MEKF_vec(dx0, P0, q_est0, omega0, nmax, F0, H, G, Q, R, t_s);

%% MMEKF Object

mmekf_obj = MM_MMEKF_vec(dx0, P0, q_est0, omega0, nmax, F0, H, G, Q, R, t_s);

%% SMEKF Object

smekf_obj = MM_SMEKF_vec(dx0, P0, q_est0, omega0, nmax, F0, H, G, Q, R, t_s);

%% Simulation
% Generate angular accelerations [x,y,z]
omega_dot = [0.000001*sin(2*pi*T./(60*60)); 0.000001*sin(2*pi*T./(60*60) - pi/2); 0.000001*cos(2*pi*T./(60*60))];

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

%% Satellite Orbit
% Costants an parameters
G = 6.674e-11;
EarthRad = 6371*1000;       % m
EarthMass = 5.97e24;        % kg
SunRad = 696340*1000;       % m

SatAltitude = 431*1000;     % m
r = EarthRad + SatAltitude; % m

% xy is the equator then this is the angle of the orbit wrt the equator.
theta = 30;                

SatVel = 2*sqrt(G*EarthMass/(r));
SatAngVel = SatVel/r;
SatPos = zeros(N,6);

% Orbit
for i = 1:N
    rad = i*t_s*SatAngVel;
    x = r*cos(rad);
    y = r*sin(rad)*cos(deg2rad(theta));
    z = -r*sin(rad)*sin(deg2rad(theta));

    rho = norm([x y z]);
    thetaE = acos(z/rho);
    phiE = 0;
    psiE = atan2(y,x);
    lat = 90-thetaE*180/pi;
    long = psiE*180/pi;
    alt = (rho - EarthRad);

    SatPos(i,:) = [x,y,z,lat,long,alt];
end

SimDate = [2022,12,9,16,20,0];

%% Star Tracker (Vector Measurements) from book example
tic
st_meas=starmeas(q_true',6,6,sig_st,nmax/nsensors)'; % (True quaternion, FOV, noise variance, max stars)
toc

%% Magnetometer
unit = symunit;
[mG2nT,~] = separateUnits(unitConvert(1*unit.mG, unit.nT));
mG2nT = double(mG2nT); % Unit conversion from milligauss to nanotesla
[nT2uT,~] = separateUnits(unitConvert(1*unit.nT, unit.uT));
nT2uT = double(nT2uT); % Unit conversion from nanotesla to millitesla
nT2mG = 1/mG2nT;
% nT2uT = nT2mG;
D1 = [0.0800, 0.0520, 0.0500;
      0.0520, 0.0500, 0.0490;
      0.0500, 0.0490, 0.0750]; % Unitless, scaling factors!
b1 = [50.00; 60.00; 55.0];
x_true = [b1; D1(1,1); D1(2,2); D1(3,3); D1(1,2); D1(1,3); D1(2,3)];
OT1 = eye(3); % ext_rot_mat_xyz(deg2rad(0.05), deg2rad(0.02), deg2rad(0.01));
Sigma1 = eye(3).*((0.15)^2);
mgm1 = Magnetometer(D1, b1, OT1, Sigma1);

% From mag_calibration_ekf_showcase
b_hat = [49.5309;
        60.8206;
         53.7601];
D_hat = [0.0810,    0.0494,    0.0495;
        0.0494,    0.0537,    0.0450;
        0.0495,    0.0450,    0.0741];

mag_meas = zeros(100,N);

tic
for i = 1:N
    [mag_true_NED] = igrfmagm(SatPos(i,6),SatPos(i,4),SatPos(i,5),decyear(SimDate)); % Magnetic field vector in North East Down frame
%     mag_true_NED = [XYZ(1);XYZ(2);-XYZ(3)]; % Magnetic field vector in North East Up frame

    rho = norm([SatPos(i,1) SatPos(i,2) SatPos(i,3)]);
    thetaE = acos(SatPos(i,3)/rho);
    psiE = atan2(SatPos(i,2),SatPos(i,1));

    mag_true_i = TIB(0,thetaE+pi,psiE)*mag_true_NED'; % Magnetic field vector in inertial frame
    
    A = quat_att_mat(q_true(:,i));
    mag_true_b = A*(mag_true_i);

    % Noisy measurements
    mag_meas1 = mgm1.simulate_reading(quat_att_mat(q_true(:,i)), mag_true_i*nT2mG); 
    mag_meas1 = (eye(3) + D_hat)*mag_meas1 - b_hat; % Calibration compensation - Should be estimated bias and scaling factors from calibration filter
    mag_meas2 = mgm1.simulate_reading(quat_att_mat(q_true(:,i)), mag_true_i*nT2mG); 
    mag_meas2 = (eye(3) + D_hat)*mag_meas2 - b_hat; % Calibration compensation
    mag_meas3 = mgm1.simulate_reading(quat_att_mat(q_true(:,i)), mag_true_i*nT2mG);
    mag_meas3 = (eye(3) + D_hat)*mag_meas3 - b_hat; % Calibration compensation

    % Normalize
    mag_meas1 = mag_meas1 / norm(mag_meas1); mag_meas2 = mag_meas2 / norm(mag_meas2); mag_meas3 = mag_meas3 / norm(mag_meas3);
    mag_true_b = (mag_true_b*nT2mG)/norm(mag_true_b*nT2mG);
    mag_true_i = (mag_true_i*nT2mG)/norm(mag_true_i*nT2mG);

    mag_meas(1:9,i) = [mag_true_b;mag_true_b;mag_true_b]; % True sun in body frame
    mag_meas(31:39,i) = [mag_meas1;mag_meas2;mag_meas3]; % Measured sun in body frame
    mag_meas(61:69,i) = [mag_true_i;mag_true_i;mag_true_i]; % True sun in inertial frame
    mag_meas(91:93,i) = 1; % Availible sensor readings
    if i == N
        display(sprintf('MAG: Done'))
    end
end
toc

%% Fine Sun Sensor
tic
global PC    % Planetary Coefficients
global const % Astronomical Constants
SAT_Const
load DE440Coeff.mat
PC = DE440Coeff;

sun_true_pos = Sun([2022,12,9,16,20,0]); % Sun in inertial frame
distE2S = norm(sun_true_pos);
sun_true_i = sun_true_pos / distE2S;
AlphaUmbra = atan2(SunRad-EarthRad,distE2S);
PV = EarthRad*distE2S/(SunRad-EarthRad);
y = sqrt(r^2 - EarthRad^2);
x = y*cos(AlphaUmbra);
UmbVert = EarthRad*(PV-x)/PV;

IsSunny = ones(1,N)*2;
sun_meas = zeros(100,N);
sig = 1*sig_sun;
for i = 1:N 
    SatVect = [SatPos(i,1) SatPos(i,2) SatPos(i,3)];
    SatHor = dot(SatVect,sun_true_pos) / dot(sun_true_pos,sun_true_pos) * sun_true_pos;
    SatHorVec = SatHor/norm(SatHor);
    ang = acos(SatHorVec(1)*sun_true_i(1) + SatHorVec(2)*sun_true_i(2) + SatHorVec(3)*sun_true_i(3));
    if(ang < 0.01)
        IsSunny(1,i) = 1;
            sun_meas(91:93,i) = 1; % Availible sensor readings
    elseif ang - pi < 0.01
        SatVert = norm(cross(SatVect,sun_true_i));
        if SatVert < UmbVert
            IsSunny(1,i) = 0;
            sun_meas(91:93,i) = 0; % Not Availible sensor readings
        else
            IsSunny(1,i) = 1;
            sun_meas(91:93,i) = 1; % Availible sensor readings
        end
    end   
    A = quat_att_mat(q_true(:,i));
    sun_true_b = A*sun_true_i; % Sun in body frame

    sun_meas1 = R_z(sig_sun*randn)*R_y(sig_sun*randn)*R_x(sig_sun*randn)*sun_true_b; sun_meas1 = sun_meas1 / norm(sun_meas1);
    sun_meas2 = R_z(sig_sun*randn)*R_y(sig_sun*randn)*R_x(sig_sun*randn)*sun_true_b; sun_meas2 = sun_meas2 / norm(sun_meas2);
    sun_meas3 = R_z(sig_sun*randn)*R_y(sig_sun*randn)*R_x(sig_sun*randn)*sun_true_b; sun_meas3 = sun_meas3 / norm(sun_meas3);

    sun_meas(1:9,i) = [sun_true_b;sun_true_b;sun_true_b]; % True sun in body frame
    sun_meas(31:39,i) = [sun_meas1;sun_meas2;sun_meas3]; % Measured sun in body frame
    sun_meas(61:69,i) = [sun_true_i;sun_true_i;sun_true_i]; % True sun in inertial frame

    if i == N
        display(sprintf('SUN: Done'))
    end
end
toc

%% Combine Sensor Measurements
meas = [
        mag_meas(1:30,:);
        sun_meas(1:30,:);
        st_meas(1:30,:);
        mag_meas(31:60,:);
        sun_meas(31:60,:);
        st_meas(31:60,:);
        mag_meas(61:90,:);
        sun_meas(61:90,:);
        st_meas(61:90,:);
        mag_meas(91:100,:);
        sun_meas(91:100,:);
        st_meas(91:100,:)];

%% Estimation time!
bias_g1 = zeros(3,N);
q_est_mekf = zeros(4, N);
q_est_norm_mekf = zeros(1,N);
bias_est_mekf = zeros(3,N);
q_est_mmekf = zeros(4, N);
q_est_norm_mmekf = zeros(1,N);
bias_est_mmekf = zeros(3,N);
q_est_smekf = zeros(4, N);
q_est_norm_smekf = zeros(1,N);
bias_est_smekf = zeros(3,N);
cond_obs_hist = zeros(1,N);
omega_meas = zeros(3,N);
%q_fu = zeros(4,N);

for i=1:N
    % Gyro measurements
    [omega_meas(:,i), bias_g1(:,i), gy1] = gy1.simulate_reading(omega(:,i));
end

tic
for i=1:N
    if i == 1;
        display(sprintf('MEKF: Started'))
    end

    if i == N/2;
        display(sprintf('MEKF: Almost there'))
    end

    % Filtering
    mekf_obj = mekf_obj.update([meas(:,i); omega_meas(:,i)]);
    mekf_obj = mekf_obj.reset();
    mekf_obj = mekf_obj.propagate();

    q_est_mekf(:,i) = mekf_obj.q_ref;
    bias_est_mekf(:,i) = mekf_obj.beta_g;
    q_est_norm_mekf(:,i) = norm(q_est_mekf(:,i));
%     cond_obs_hist(:,i) = mekf_obj.cond_obs;
        if i == N;
        display(sprintf('MEKF: Finished'))
    end
end
toc

tic
for i=1:N
    if i == 1;
        display(sprintf('MMEKF: Started'))
    end
    % Gyro measurements

    if i == N/2;
        display(sprintf('MMEKF: Almost there'))
    end

    % Filtering
    mmekf_obj = mmekf_obj.update([meas(:,i); omega_meas(:,i)]);
    mmekf_obj = mmekf_obj.reset();
    mmekf_obj = mmekf_obj.propagate();

    q_est_mmekf(:,i) = mmekf_obj.q_ref;
    bias_est_mmekf(:,i) = mmekf_obj.beta_g;
    q_est_norm_mmekf(:,i) = norm(q_est_mmekf(:,i));
%     cond_obs_hist(:,i) = mekf_obj.cond_obs;
        if i == N;
        display(sprintf('MMEKF: Finished'))
    end
end
toc

tic
for i=1:N
    if i == 1;
        display(sprintf('SMEKF: Started'))
    end

    if i == N/2;
        display(sprintf('SMEKF: Almost there'))
    end

    % Filtering
    smekf_obj = smekf_obj.update([meas(:,i); omega_meas(:,i)]);
    smekf_obj = smekf_obj.reset();
    smekf_obj = smekf_obj.propagate();

    q_est_smekf(:,i) = smekf_obj.q_ref;
    bias_est_smekf(:,i) = smekf_obj.beta_g;
    q_est_norm_smekf(:,i) = norm(q_est_smekf(:,i));
%     cond_obs_hist(:,i) = mekf_obj.cond_obs;
        if i == N;
        display(sprintf('SMEKF: Finished'))
    end
end
toc

% Compute estimation errors
err_mekf = zeros(3, N);
for i=1:N
    err_mekf(:,i) = quat2eul(quaternion(q_true(:,i)'))' - quat2eul(quaternion(q_est_mekf(:,i)'))';
end

err_mmekf = zeros(3, N);
for i=1:N
    err_mmekf(:,i) = quat2eul(quaternion(q_true(:,i)'))' - quat2eul(quaternion(q_est_mmekf(:,i)'))';
end

err_smkef = zeros(3, N);
for i=1:N
    err_smekf(:,i) = quat2eul(quaternion(q_true(:,i)'))' - quat2eul(quaternion(q_est_smekf(:,i)'))';
end

% Compute estimation errors
ag_err_mekf = zeros(3, N);
for i=1:N
    ag_err_mekf(:,i) = quat2gvr(sh_quat_mult(q_true(:,i), quat_inv(q_est_mekf(:,i))));
end
err_mekf = ag_err_mekf;

ag_err_mmekf = zeros(3, N);
for i=1:N
    ag_err_mmekf(:,i) = quat2gvr(sh_quat_mult(q_true(:,i), quat_inv(q_est_mmekf(:,i))));
end
err_mmekf = ag_err_mmekf;

ag_err_smekf = zeros(3, N);
for i=1:N
    ag_err_smekf(:,i) = quat2gvr(sh_quat_mult(q_true(:,i), quat_inv(q_est_smekf(:,i))));
end
err_smekf = ag_err_smekf;

%% Plot-time!
T_plot = T*t_s;
PLT_CFG = [0; 0; 1; 1; 1; 0; 0; 0];
bias = bias_g1;

% Plot angular accelerations
if PLT_CFG(1)
figure;
plot(T_plot, omega_dot, "LineWidth", 1)
legend('$\dot{\omega}_x$', '$\dot{\omega}_y$', '$\dot{\omega}_z$', "Interpreter", "latex", "FontSize", 14);
end

% Plot resulting angular velocities
if PLT_CFG(2)
figure;
plot(T_plot, omega)
legend('$\omega_x$', '$\omega_y$', '$\omega_z$', "Interpreter", "latex", "FontSize", 14);
end

% % Plot estimation errors in Euler angles
% if PLT_CFG(3)
% figure; hold on; sgtitle('MEKF: Estimation errors in Euler angles')
% plot(T_plot, err);
% legend('$\epsilon_z$', '$\epsilon_y$', '$\epsilon_x$', "Interpreter", "latex", "FontSize", 14);
% end

% % Plot true attitude over estimated attitude
% if PLT_CFG(4)
% figure; hold on; sgtitle('MEKF: True attitude vs estimated attitude')
% subplot(2,2,1)
% plot(T_plot, q_est(1,:), T_plot, q_true(1,:));
% legend('$\hat{q_1}$', '$q_1$', "Interpreter", "latex", "FontSize", 10);
% subplot(2,2,2)
% plot(T_plot, q_est(2,:), T_plot, q_true(2,:));
% legend('$\hat{q_2}$', '$q_2$', "Interpreter", "latex", "FontSize", 10);
% subplot(2,2,3)
% plot(T_plot, q_est(3,:), T_plot, q_true(3,:));
% legend('$\hat{q_3}$', '$q_3$', "Interpreter", "latex", "FontSize", 10);
% subplot(2,2,4)
% plot(T_plot, q_est(4,:), T_plot, q_true(4,:));
% legend('$\hat{q_4}$', '$q_4$', "Interpreter", "latex", "FontSize", 10);
% end

% Plot biases
% if PLT_CFG(5)
% figure; hold on; sgtitle('MEKF: Biases')
% subplot(3,1,1)
% plot(T_plot, bias(1,:)*180/pi*3600, T_plot, bias_est(1,:)*180/pi*3600)
% xlim([T_plot(:,1), T_plot(:,end)])
% % ylim([-0.1, 0.2])
% legend('${\beta}_x$', '$\hat{\beta}_x$', "Interpreter", "latex", "FontSize", 14);
% subplot(3,1,2)
% plot(T_plot, bias(2,:)*180/pi*3600, T_plot, bias_est(2,:)*180/pi*3600)
% xlim([T_plot(1), T_plot(end)])
% legend('${\beta}_y$', '$\hat{\beta}_y$', "Interpreter", "latex", "FontSize", 14);
% % ylim([-0.1, 0.2])
% subplot(3,1,3)
% plot(T_plot, bias(3,:)*180/pi*3600, T_plot, bias_est(3,:)*180/pi*3600)
% xlim([T_plot(1), T_plot(end)])
% % ylim([-0.1, 0.2])
% legend('${\beta}_z$', '$\hat{\beta}_z$', "Interpreter", "latex", "FontSize", 14);
% end

% Plot measured vs true velocities
if PLT_CFG(6)
figure; hold on;
subplot(3,1,1)
plot(T_plot, omega_meas(1,:)-omega(1,:))
subplot(3,1,2)
plot(T_plot, omega(2,:)-omega_meas(2,:))
subplot(3,1,3)
plot(T_plot, omega(3,:)-omega_meas(3,:))
end

if PLT_CFG(7)
figure; hold on; sgtitle('Norm and conditioning number')
subplot(2,1,1)
plot(T_plot, q_norm, T_plot, q_est_norm)
subplot(2,1,2)
plot(T_plot,cond_obs_hist)
end


% Plot estimation errors in Euler angles
if PLT_CFG(3)
figure; hold on; sgtitle('Estimation errors in Euler angles')
plot(T_plot, err_mekf(:,:),'r',T_plot, err_mmekf(:,:),'g',T_plot, err_smekf(:,:),'b');
legend('$MEKF error$','$MMEKF error$','$SMEKF error$',"Interpreter", "latex", "FontSize", 10);
end

% Plot true attitude over estimated attitude
if PLT_CFG(4)
figure; hold on; sgtitle('MEKF: True attitude vs estimated attitude')
subplot(2,2,1)
plot(T_plot, q_true(1,:), T_plot, q_est_mekf(1,:), T_plot, q_est_mmekf(1,:), T_plot, q_est_smekf(1,:));
legend('$q_1$','$\hat{q_1}_{mekf}$','$\hat{q_1}_{mmekf}$', '$\hat{q_1}_{smekf}$', "Interpreter", "latex", "FontSize", 10);
subplot(2,2,2)
plot(T_plot, q_true(2,:), T_plot, q_est_mekf(2,:), T_plot, q_est_mmekf(2,:), T_plot, q_est_smekf(2,:));
legend('$q_2$','$\hat{q_2}_{mekf}$','$\hat{q_2}_{mmekf}$','$\hat{q_2}_{smekf}$',"Interpreter", "latex", "FontSize", 10);
subplot(2,2,3)
plot(T_plot, q_true(3,:), T_plot, q_est_mekf(3,:),T_plot, q_est_mmekf(3,:),T_plot, q_est_smekf(3,:));
legend('$q_3$', '$\hat{q_3}_{mekf}$', '$\hat{q_3}_{mmekf}$', '$\hat{q_3}_{smekf}$', "Interpreter", "latex", "FontSize", 10);
subplot(2,2,4)
plot(T_plot, q_true(4,:),T_plot, q_est_mekf(4,:),T_plot, q_est_mmekf(4,:),T_plot, q_est_smekf(4,:));
legend('$q_4$','$\hat{q_4}_{mekf}$','$\hat{q_4}_{mmekf}$','$\hat{q_4}_{smekf}$', "Interpreter", "latex", "FontSize", 10);
end

if PLT_CFG(5)
figure; hold on; sgtitle('Biases')
subplot(3,1,1)
plot(T_plot, bias(1,:)*180/pi*3600, T_plot, bias_est_mekf(1,:)*180/pi*3600,T_plot, bias_est_mmekf(1,:)*180/pi*3600,T_plot, bias_est_smekf(1,:)*180/pi*3600)
xlim([T_plot(:,1), T_plot(:,end)])
% ylim([-0.1, 0.2])
legend('${\beta}_x$', '$\hat{\beta}_x^{MEKF}$','$\hat{\beta}_x^{MMEKF}$','$\hat{\beta}_x^{SMEKF}$', "Interpreter", "latex", "FontSize", 14);
subplot(3,1,2)
plot(T_plot, bias(2,:)*180/pi*3600, T_plot, bias_est_mekf(2,:)*180/pi*3600,T_plot, bias_est_mmekf(2,:)*180/pi*3600,T_plot, bias_est_smekf(2,:)*180/pi*3600)
xlim([T_plot(1), T_plot(end)])
legend('${\beta}_y$', '$\hat{\beta}_y^{MEKF}$','$\hat{\beta}_y^{MMEKF}$','$\hat{\beta}_y^{SMEKF}$', "Interpreter", "latex", "FontSize", 14);
% ylim([-0.1, 0.2])
subplot(3,1,3)
plot(T_plot, bias(3,:)*180/pi*3600, T_plot, bias_est_mekf(3,:)*180/pi*3600, T_plot, bias_est_mmekf(3,:)*180/pi*3600, T_plot, bias_est_smekf(3,:)*180/pi*3600)
xlim([T_plot(1), T_plot(end)])
% ylim([-0.1, 0.2])
legend('${\beta}_z$', '$\hat{\beta}_z^{MEKF}$', '$\hat{\beta}_z^{MMEKF}$', '$\hat{\beta}_z^{SMEKF}$', "Interpreter", "latex", "FontSize", 14);
end


mse_mekf = [mean(err_mekf(1,1:end-1).^2), mean(err_mekf(2,1:end-1).^2), mean(err_mekf(3,1:end-1).^2)]
bias_err_mekf = (bias(:,1000:N) - bias_est_mekf(:,1000:N))*180/pi*3600;
bias_mse_mekf = mean(bias_err_mekf.^2, 2)

mse_mmekf = [mean(err_mmekf(1,1:end-1).^2), mean(err_mmekf(2,1:end-1).^2), mean(err_mmekf(3,1:end-1).^2)]
bias_err_mmekf = (bias(:,1000:N) - bias_est_mmekf(:,1000:N))*180/pi*3600;
bias_mse_mmekf = mean(bias_err_mmekf.^2, 2)

mse_smekf = [mean(err_smekf(1,1:end-1).^2), mean(err_smekf(2,1:end-1).^2), mean(err_smekf(3,1:end-1).^2)]
bias_err_smekf = (bias(:,1000:N) - bias_est_smekf(:,1000:N))*180/pi*3600;
bias_mse_smekf = mean(bias_err_smekf.^2, 2)









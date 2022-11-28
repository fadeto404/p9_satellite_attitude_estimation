close all; clear; %clc;
%% MEKF showcase simulation with fused star tracker measurement and gyro
rng(42);
% Simulation time basics
t_s = 0.1; % [s]
num_min = 90; % Sim time in [min]
N = num_min*60 / t_s; % No. of time steps
T = 0:t_s:(N-1)*t_s;

%% Star tracker setup
% Rotation order: Z-Y-X
variances = [deg2rad(360/(60^2)), deg2rad(36/(60^2)), deg2rad(36/(60^2))];
variances2 = [deg2rad(36/(60^2)), deg2rad(360/(60^2)), deg2rad(36/(60^2))];
R1 = diag(variances); % Euler angle covariance matrix
R2 = diag(variances2);
% BigR(:,:,1) = R1;
% BigR(:,:,2) = R2;
% 
% st1 = StarTracker(variances);
% st2 = StarTracker(variances2);
% quat_fuser = QuaternionFuser(BigR);

%% New star tracker setup (preferred)
% Variance of angular error around X, Y, and Z in the local star tracker
% frame
variances = [deg2rad(36/(60^2)), deg2rad(36/(60^2)), deg2rad(360/(60^2))];
R_st = diag(variances); % Star tracker frame noise covariance matrix
ST1_rot = eye(3); % Star tracker 1 rotation wrt. body frame
ST2_rot = [1, 0, 0; 0, 0, -1; 0, 1, 0]; % Star tracker 2 rotation wrt. BF
st1 = StarTrackerSimple(variances, ST1_rot); % Star tracker 1 init
st2 = StarTrackerSimple(variances, ST2_rot); % Star tracker 2 init

BigR(:,:,1) = ST1_rot*R_st; % ST1 covariance in space craft body frame
BigR(:,:,2) = ST2_rot*R_st; % ST2 ...
quat_fuser = QuaternionFuser(BigR); % Quaternion averaging init

%% Gyro setup
b0 = ones(3,1)*(0.1 / (180/pi*3600)); % 0.1 [deg/hr] to [rad/s] (init bias)
gyro_noise_variance = (sqrt(10)*10^(-7))^2; % Noise variance
gyro_bias_variance = (sqrt(10)*10^(-10))^2; % Bias variance
gyro_noise_cov = eye(3)*gyro_noise_variance; % Noise covariance matrix
gyro_bias_cov = eye(3)*gyro_bias_variance; % Bias random walk cov. mat.
gy1 = Gyro(gyro_noise_variance, gyro_bias_variance, b0, t_s); % Gyro sim init

%% Filter setup
dx0 = [zeros(3,1); zeros(3,1)];
P0 = [(0.1*pi/180)^2*eye(3), zeros(3); zeros(3), (0.2*pi/180/3600)^2*eye(3)];
P0 = blkdiag((0.1*pi/180/3600)^2*eye(3) * eye(3), (0.2*pi/180/3600)^2 * eye(3));
q_est0 = 1/sqrt(2)*[1; 0; 0; 1];
omega0 = 0.0000001*ones(3,1);
% omega0 = [inv(91.5/(2*pi)*60)*0; -inv(91.5/(2*pi)*60); 0]; % For constant velocity sim
% omega0 = 0.1*pi/180*[sin(t_s*0.01*T(1)') sin(t_s*0.0085*T(1)') cos(t_s*0.0085*T(1)')]';
F0 = [zeros(3) -eye(3); zeros(3,6)];
G = [-eye(3) zeros(3); zeros(3) eye(3)];
H = [eye(3) zeros(3)];
% R = quat_fuser.R_bar;
R = R_st;
% R = diag([6/3600*pi/180, 6/3600*pi/180, 6/3600*pi/180]).^2;
Q = [gyro_noise_cov, zeros(3); zeros(3), gyro_bias_cov];

mekf_obj = MM_MEKF(dx0, P0, q_est0, omega0, F0, H, G, Q, R, t_s);

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

% Estimation time!
q_est = zeros(4, N);
q_est_norm = zeros(1,N);
bias_g1 = zeros(3,N);
bias_est = zeros(3,N);
cond_obs_hist = zeros(1,N);
omega_meas = zeros(3,N);
%q_fu = zeros(4,N);

for i=1:N
    % Gyro measurements
    [omega_meas(:,i), bias_g1(:,i), gy1] = gy1.simulate_reading(omega(:,i));

    % Quaternion measurements
    q1 = st1.simulate_reading(q_true(:,i));
    q2 = st2.simulate_reading(q_true(:,i));
    quat_fuser = quat_fuser.fuse([q1, q2]);
    [q1, q1err] = st1.simulate_reading(q_true(:,i));
    [q2, q2err] = st2.simulate_reading(q_true(:,i));
    %quat_fuser = quat_fuser.fuse([q1, q2]);
    %q_bar = quat_fuser.q_bar;
    q_bar = q1;
    
    % Filtering
    mekf_obj = mekf_obj.update([q_bar; omega_meas(:,i)]);
    mekf_obj = mekf_obj.reset();
    mekf_obj = mekf_obj.propagate();
    q_est(:,i) = mekf_obj.q_ref;
    bias_est(:,i) = mekf_obj.beta_g;
    q_est_norm(:,i) = norm(q_est(:,i));
    %cond_obs_hist(:,i) = mekf_obj.cond_obs;
end

% Compute estimation errors
err = zeros(3, N);
for i=1:N
    err(:,i) = quat2eul(quaternion(q_true(:,i)'))' - quat2eul(quaternion(q_est(:,i)'))';
end

% Compute estimation errors
ag_err = zeros(3, N);
for i=1:N
    ag_err(:,i) = quat2gvr(sh_quat_mult(q_true(:,i), quat_inv(q_est(:,i))));
end
err = ag_err;

%% Plot-time!
T_plot = T*t_s;
PLT_CFG = [0; 1; 1; 1; 1; 1; 0];
% Plot angular accelerations
if PLT_CFG(1)
figure(1);
plot(T_plot, omega_dot, "LineWidth", 1)
legend('$\dot{\omega}_x$', '$\dot{\omega}_y$', '$\dot{\omega}_z$', "Interpreter", "latex", "FontSize", 14);
end

% Plot resulting angular velocities
if PLT_CFG(2)
figure(2);
plot(T_plot, omega)
legend('$\omega_x$', '$\omega_y$', '$\omega_z$', "Interpreter", "latex", "FontSize", 14);
end

% Plot estimation errors in Euler angles
if PLT_CFG(3)
figure(3); hold on;
plot(T_plot, err);
legend('$\epsilon_z$', '$\epsilon_y$', '$\epsilon_x$', "Interpreter", "latex", "FontSize", 14);
end

% Plot true attitude over estimated attitude
if PLT_CFG(4)
figure(4); hold on;
subplot(2,2,1)
plot(T_plot, q_est(1,:), T_plot, q_true(1,:));
legend('$\hat{q_1}$', '$q_1$', "Interpreter", "latex", "FontSize", 10);
subplot(2,2,2)
plot(T_plot, q_est(2,:), T_plot, q_true(2,:));
legend('$\hat{q_2}$', '$q_2$', "Interpreter", "latex", "FontSize", 10);
subplot(2,2,3)
plot(T_plot, q_est(3,:), T_plot, q_true(3,:));
legend('$\hat{q_3}$', '$q_3$', "Interpreter", "latex", "FontSize", 10);
subplot(2,2,4)
plot(T_plot, q_est(4,:), T_plot, q_true(4,:));
legend('$\hat{q_4}$', '$q_4$', "Interpreter", "latex", "FontSize", 10);
end

% Plot biases
if PLT_CFG(5)
bias = bias_g1;
figure(5); hold on;
subplot(3,1,1)
plot(T_plot, bias(1,:)*180/pi*3600, T_plot, bias_est(1,:)*180/pi*3600)
xlim([T_plot(:,1), T_plot(:,end)])
% ylim([-0.1, 0.2])
legend('${\beta}_x$', '$\hat{\beta}_x$', "Interpreter", "latex", "FontSize", 14);
subplot(3,1,2)
plot(T_plot, bias(2,:)*180/pi*3600, T_plot, bias_est(2,:)*180/pi*3600)
xlim([T_plot(1), T_plot(end)])
legend('${\beta}_y$', '$\hat{\beta}_y$', "Interpreter", "latex", "FontSize", 14);
% ylim([-0.1, 0.2])
subplot(3,1,3)
plot(T_plot, bias(3,:)*180/pi*3600, T_plot, bias_est(3,:)*180/pi*3600)
xlim([T_plot(1), T_plot(end)])
% ylim([-0.1, 0.2])
legend('${\beta}_z$', '$\hat{\beta}_z$', "Interpreter", "latex", "FontSize", 14);
end

% Plot measured vs true velocities
if PLT_CFG(6)
figure(6); hold on;
subplot(3,1,1)
plot(T_plot, omega_meas(1,:)-omega(1,:))
subplot(3,1,2)
plot(T_plot, omega(2,:)-omega_meas(2,:))
subplot(3,1,3)
plot(T_plot, omega(3,:)-omega_meas(3,:))
end

if PLT_CFG(7)
figure(7); hold on;
subplot(2,3,4)
plot(T_plot, q_norm, T_plot, q_est_norm)
subplot(2,3,5)
plot(T_plot,cond_obs_hist)
end

mse = [mean(err(1,1:end-1).^2), mean(err(2,1:end-1).^2), mean(err(3,1:end-1).^2)]
bias_err = (bias - bias_est)*180/pi*3600;
bias_mse = mean(bias_err.^2, 2)

















clear; clc; close all;
%% MEKF showcase simulation with fused star tracker measurement and gyro
% Simulation time basics
t_s = 0.01; % [s]
num_min = 60; % Sim time in [min]
N = num_min*60 / t_s; % No. of time steps
T = 0:t_s:(N-1)*t_s;

%% Star tracker setup
% Rotation order: Z-Y-X
variances = [deg2rad(50/(60^2)), deg2rad(5/(60^2)), deg2rad(5/(60^2))];
variances2 = [deg2rad(5/(60^2)), deg2rad(50/(60^2)), deg2rad(5/(60^2))];
R1 = diag(variances); % Euler angle covariance matrix
R2 = diag(variances2);
BigR(:,:,1) = R1;
BigR(:,:,2) = R2;

st1 = StarTracker(variances);
st2 = StarTracker(variances2);
quat_fuser = QuaternionFuser(BigR);

%% Gyro setup
gyro_noise_cov = eye(3)*(sqrt(10)*10^(-7))^2; % TODO: Find reasonable cov matrices
gyro_bias_cov = eye(3)*(sqrt(10)*10^(-10))^2;
gy1 = Gyro(gyro_noise_cov, gyro_bias_cov, zeros(3,1), t_s);

%% Filter setup
x0 = zeros(6,1);
P0 = [(0.1*pi/180)^2*eye(3), zeros(3); zeros(3), (0.2*pi/180/3600)^2*eye(3)];
q_est0 = 1/sqrt(2)*[1; 0; 0; 1];
omega0 = 0.0000001*ones(3,1);
F0 = [zeros(3) -eye(3); zeros(3,6)];
G = [-eye(3) zeros(3); zeros(3) eye(3)];
H = [eye(3) zeros(3)];
R = quat_fuser.R_bar*100;
Q = [gyro_noise_cov, zeros(3); zeros(3), gyro_bias_cov];
mekf_obj = MM_MEKF(x0, P0, q_est0, omega0, F0, H, G, Q, R, t_s);

%% Simulation
% Generate angular accelerations [x,y,z]
omega_dot = [0.0001*sin(2*pi*T./(60*15)); 0.0001*sin(2*pi*T./(60*15) - pi/2); 0.0001*cos(2*pi*T./(60*15))];

% Generate angular velocities by Euler integration
omega = zeros(3,N);
for i=1:N-1
    omega(:,i+1) = omega(:,i) + t_s*omega_dot(:,i);
end

% Compute true attitude quaternion
q0 = 1/sqrt(2)*[1; 0; 0; 1];
%q0 = randn(4,1);
%q0 = q0./norm(q0);
q_true = [q0, zeros(4,N-1)];
q_norm = [1, zeros(1,N-1)];

% Integrating quaternions
% Convert ang. vel. to unit quaternion, multiply prev. q_true by q_omega
for i=1:N-1
    %A = 0.5*sh_prod_mat([omega(:,i); 0]);
    if omega(:,i) == zeros(3,1)
        q_true(:,i+1) = q_true(:,i);
        q_norm(:,i+1) = norm(q_true(:,i+1));
    else
        omega_norm = norm(omega(:,i));
        omega_normalized = omega(:,i)./omega_norm;
        q_omega = [omega_normalized*sin(0.5*omega_norm*t_s); cos(0.5*omega_norm*t_s)];
        q_true(:,i+1) = sh_quat_mult(q_omega, q_true(:,i));
        q_norm(:,i+1) = norm(q_true(:,i+1));
    end
end

% Estimation time!
q_est = zeros(4, N);
q_est_norm = zeros(1,N);
bias = zeros(3,N);
bias_est = zeros(3,N);
for i=1:N-1
    % Gyro measurements
    [omega_meas, bias(:,i), gy1] = gy1.simulate_reading(omega(:,i));

    % Quaternion measurements
    q1 = st1.simulate_reading(quaternion(q_true(:,i)'));
    q2 = st2.simulate_reading(quaternion(q_true(:,i)'));
    quat_fuser = quat_fuser.fuse([q1.compact', q2.compact']);
    q_bar = quat_fuser.q_bar;

    mekf_obj = mekf_obj.propagate();
    mekf_obj = mekf_obj.update([q_bar; omega_meas]);
    mekf_obj = mekf_obj.reset();
    q_est(:,i) = mekf_obj.q_ref;
    bias_est(:,i) = mekf_obj.beta_g;
    q_est_norm(:,i) = norm(q_est(:,i));
end

% Compute estimation errors
err = zeros(3, N);
for i=1:N
    err(:,i) = quat2eul(quaternion(q_true(:,i)'))' - quat2eul(quaternion(q_est(:,i)'))';
end

%% Plot-time!
% Plot angular accelerations
figure(1);
plot(T, omega_dot, "LineWidth", 1)
legend('$\dot{\omega}_x$', '$\dot{\omega}_y$', '$\dot{\omega}_z$', "Interpreter", "latex", "FontSize", 14);

% Plot resulting angular velocities
figure(2);
plot(T, omega)
legend('$\omega_x$', '$\omega_y$', '$\omega_z$', "Interpreter", "latex", "FontSize", 14);

% Plot estimation errors in Euler angles
figure(3); hold on;
plot(T, err);
legend('$\epsilon_z$', '$\epsilon_y$', '$\epsilon_x$', "Interpreter", "latex", "FontSize", 14);

% Plot true attitude over estimated attitude
figure(4); hold on;
subplot(2,2,1)
plot(T, q_est(1,:), T, q_true(1,:));
legend('$\hat{q_1}$', '$q_1$', "Interpreter", "latex", "FontSize", 10);
subplot(2,2,2)
plot(T, q_est(2,:), T, q_true(2,:));
legend('$\hat{q_2}$', '$q_2$', "Interpreter", "latex", "FontSize", 10);
subplot(2,2,3)
plot(T, q_est(3,:), T, q_true(3,:));
legend('$\hat{q_3}$', '$q_3$', "Interpreter", "latex", "FontSize", 10);
subplot(2,2,4)
plot(T, q_est(4,:), T, q_true(4,:));
legend('$\hat{q_4}$', '$q_4$', "Interpreter", "latex", "FontSize", 10);

% Plot various other info
figure(5); hold on;
subplot(3,3,1)
plot(T, bias(1,:)*180/pi*3600, T, bias_est(1,:)*180/pi*3600)
subplot(3,3,2)
plot(T, bias(2,:)*180/pi*3600, T, bias_est(2,:)*180/pi*3600)
subplot(3,3,3)
plot(T, bias(3,:)*180/pi*3600, T, bias_est(3,:)*180/pi*3600)
subplot(3,3,4)
plot(T, q_norm, T, q_est_norm)


mse = [mean(err(1,1:end-1).^2), mean(err(2,1:end-1).^2), mean(err(3,1:end-1).^2)]

















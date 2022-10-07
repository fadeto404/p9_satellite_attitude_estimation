%% Kinematics simulation with constant angular velocity of PCS in LVLH
% Quaternions
q0 = quaternion(1,0,0,0);
%q1 = quaternion(0.5, 0.7, 0.8, 0.5)./norm(quaternion(0.5, 0.7, 0.8, 0.5));
%q2 = quaternion(0.3, 0.2, 0.9, 0.8)./norm(quaternion(0.3, 0.2, 0.9, 0.8));

Omega_po = eul2quat([0 1 0])
Omega_po = quaternion(Omega_po(1), Omega_po(2), Omega_po(3), Omega_po(4)).normalize()
%Omega_po = quaternion(0, 1, 0.7, 0.5).normalize();
Omega_po_eul = quat2eul(Omega_po)
[omega0, omega1, omega2, omega3] = Omega_po.parts;
%Delta_phi = quaternion(0.1,0,0,0);
%epsilon = quaternion(0,0,1,0);
R = @(Omega)([0 -omega1 -omega2 -omega3; ...
              omega1 0 -omega3 omega2; ...
              omega2 omega3 0 -omega1; ...
              omega3 -omega2 omega1 0]);
q_dot_po = @(q, t)(1/2 * R(Omega_po) * q);



% Euler integration
N = 10000 % 10000 ms = 10s
t_s = 1/1000; % 1 ms
q_dot_sim = [];
q_sim = [q0.compact];
for i = 1:N
    q_dot_sim(i,:) = q_dot_po(q_sim(i,:)', i); % 1/2 * R(Omega_po) * (q_sim(i,:))';
    q_sim(i+1,:) = q_sim(i,:) + t_s*q_dot_sim(i,:);
end

figure(1);
subplot(2,1,1);
plot(0:N,q_sim);
legend('q0', 'q1', 'q2', 'q3')
xlim([0, N]);
subplot(2,1,2);
plot(0:N, quat2eul(q_sim))

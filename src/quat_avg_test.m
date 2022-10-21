clear; clc; close all;
% Code showcasing quaternion fusion by averaging, weighted averaging, or 
% selecting
rng(69420);

%% Generate test quaternions
q_true = randn(4,1);
q_true = q_true./norm(q_true);

% Rotation order: Z-Y-X
variances = [deg2rad(50/(60^2)), deg2rad(5/(60^2)), deg2rad(5/(60^2))];
variances2 = [deg2rad(5/(60^2)), deg2rad(50/(60^2)), deg2rad(5/(60^2))];
R1 = diag(variances); % Euler angle covariance matrix
R2 = diag(variances2);
BigR(:,:,1) = R1; BigR(:,:,2) = R2;

st1 = StarTracker(variances);
st2 = StarTracker(variances2);
QF = QuaternionFuser(BigR);
for i=1:25000
    q1 = st1.simulate_reading(quaternion(q_true'));
    q2 = st2.simulate_reading(quaternion(q_true'));
    q_avg_matlab = meanrot([q1;q2]).compact';
    q1 = q1.compact';
    q2 = q2.compact';
    
    %% Quaternion averaging; scalar weights
    w = [0.5; 0.5]; % w1 = w2 ---> Simple average, equivalent to meanrot() fct.
    w1 = w(1);
    w2 = w(2);
    z = sqrt((w(1)-w(2))^2 + 4*w(1)*w(2)*(q1'*q2)^2);
    q_avg = sqrt((w1* (w1-w2+z))/(z*(w1+w2+z)))*q1 + sign(q1'*q2)*sqrt((w2* (w2-w1+z))/(z*(w1+w2+z)))*q2;
    
    %% Quaternion averaging; matrix weights using covariance matrices -> MLE
    %R1 = diag(variances); % Euler angle covariance matrix
    %R2 = diag(variances2);
    
    %R1 = 2*eye(3); % Equivalent to the simple average
    %R2 = R1;
    
    M = -(quat_xi_mat(q1)*inv(R1)*(quat_xi_mat(q1)') + quat_xi_mat(q2)*inv(R2)*(quat_xi_mat(q2)'));
    [eigvec, eigval] = eig(M);
    [max_eigval, ind] = max(diag(eigval));
    q_bar = eigvec(:,ind);
    
    % Same process, just by using the class instead
    QF = QF.fuse([q1, q2]);
    q_bar2 = QF.q_bar;
    
    %% Quaternion "average" by selecting low-noise axes
    % Downside: Requires converting to Euler angles to pick out the axes
    eul1 = quat2eul(q1', "ZYX")';
    eul2 = quat2eul(q2', "ZYX")';
    eul_select = [eul2(1); eul1(2); eul1(3)];
    q_select = eul2quat(eul_select', "ZYX")';
    
    %% Error signal computation
    errq1 = sh_quat_mult(q_true,quat_inv(q_bar)); errq1([1 2 3 4]) = errq1([4 1 2 3]);
    errq2 = sh_quat_mult(q_true,quat_inv(q_avg)); errq2([1 2 3 4]) = errq2([4 1 2 3]);
    errq3 = sh_quat_mult(q_true,quat_inv(q_avg_matlab)); errq3([1 2 3 4]) = errq3([4 1 2 3]);
    errq4 = sh_quat_mult(q_true,quat_inv(q_select)); errq4([1 2 3 4]) = errq4([4 1 2 3]);
    errq5 = sh_quat_mult(q_true,quat_inv(q1)); errq5([1 2 3 4]) = errq5([4 1 2 3]); % If we just used q1, ignoring star tracker 2
    errq6 = sh_quat_mult(q_true,quat_inv(q_bar2)); errq6([1 2 3 4]) = errq6([4 1 2 3]);
    err1(i,:) = quat2eul(errq1');
    err2(i,:) = quat2eul(errq2');
    err3(i,:) = quat2eul(errq3');
    err4(i,:) = quat2eul(errq4');
    err5(i,:) = quat2eul(errq5');
    err6(i,:) = quat2eul(errq6');

end

%% Evaluation/comparison (MSE for each Euler angle axis)
mle = (mean(err1.^2, 1))
avg = (mean(err2.^2, 1))
avg2 = (mean(err3.^2, 1))
selection = (mean(err4.^2, 1))
st1_pref = (mean(err5.^2, 1))
QF_err = (mean(err6.^2, 1))


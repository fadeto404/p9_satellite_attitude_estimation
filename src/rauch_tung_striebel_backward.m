% Rauch-Tung-Striebel Kalman smoother backward pass
% All arguments given must be a batch of saved vectors/matrices from the
% forward pass (the EKF) of N samples. States are assumed rows, time
% columns in state vector x.
function [x_s, P_s] = rauch_tung_striebel_backward(x_pre, x_post, P_pre, P_post, Phi)
    [n, N] = size(x_post); % Num. states, num. samples
    x_s = zeros(n, N);      % Allocate size
    P_s = zeros(n, n, N);   % Allocate size
    
    x_s(:,N) = x_post(:,N); % Initialize smoothed state (starts at back)
    P_s(:,:,N) = P_post(:,:,N); % Initialize smoothed error covariance
    for i=N-1:-1:1
        K_k = P_post(:,:,i) * Phi(:,:,i)' / (P_pre(:,:,i+1)); % Smoothing gain
        x_s(:,i) = x_post(:,i) + K_k*(x_s(:,i+1) - x_pre(:,i+1)); % State smoothing
        P_s(:,:,i) = P_post(:,:,i) - K_k*(P_pre(:,:,i+1) - P_s(:,:,i+1))*K_k'; % Covariance smoothing
    end
end
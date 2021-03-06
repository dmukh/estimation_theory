% ESTIMATION FINAL PROJECT - EXTENDED SEQUENTIAL
clear;

% Read data
load obs.mat
tVec = obs(:,1);
rangeVec = obs(:,2);
weightVec = obs(:,3);
sIDvec = obs(:,4);

% Constants
global re h0 gamma0 rho0 omegaE theta0
re = 6378136.3;
h0 = 9.2e5;
gamma0 = 5.381e-6;
rho0 = 4.36e-14;
omegaE = 7.2921157746e-5;
theta0 = 1.6331958133;

% Initial values
X0=[4973900 -4300600 3486200 2850 5820 3470 3.986004415e14 1.1926268e-3 .0375 -2517400 -4198500 4076500];
xbar = zeros(12, 1);
sigmas = [100 100 100 10 10 10 1e6 1e-2 .01 100 100 100];
Pbar0 = diag(sigmas.^2);
phi0 = eye(12);

options = odeset('RelTol',1e-5);

tic
% Start filter
ssrPre = 0;
ssrPost = 0;
stationDat = cell(3,1);
stationDatPost = cell(3,1);
P = Pbar0;
xhat = xbar;

for i = 1:length(tVec)
    if i == 1
        X_i = X0;
        phi_i = phi0;
    elseif tVec(i) == tVec(i-1)
        phi_i = phi0;
    elseif tVec(i) ~= tVec(i-1)
        % Propagate next time span
        X0comb = zeros(156, 1);
        X0comb(1:12) = X_i;
        X0comb(13:end) = reshape(phi0, 144, 1);
        [~, Xout] = ode45(@populateXdot, [tVec(i-1) tVec(i)], X0comb, options);
        X_i = reshape(Xout(end,1:12), 12, 1)';
        phi_i = reshape(Xout(end,13:end), 12, 12);
    end
    
    % Time update
    Pbar = phi_i * P * phi_i';
    
    % Obs and gain
    sID = sIDvec(i);
    [Htilde, G] = populateMeas(tVec(i), X_i);
    yi = rangeVec(i) - G(sID);
    Ki = Pbar*Htilde(sID,:)' * inv(Htilde(sID,:)*Pbar*Htilde(sID,:)'+1/weightVec(i));
    
    % Measurement and orbit update
    xhat = Ki * yi;
    X_i = X_i + xhat';
    P = (eye(12) - Ki*Htilde(sID,:)) * Pbar;
    
    stationDat{sID}(end+1,:) = [tVec(i) yi];
    ssrPre = ssrPre + yi^2;
    
    % Post fit
    [~, Gpost] = populateMeas(tVec(i), X_i);
    yiPost = rangeVec(i) - Gpost(sID);
    ssrPost = ssrPost + yiPost^2;
    stationDatPost{sID}(end+1,:) = [tVec(i) yiPost];
    
end

% Backwards mapping
X0comb = zeros(156, 1);
X0comb(1:12) = X_i;
X0comb(13:end) = reshape(phi0, 144, 1);
[~, Xout] = ode45(@populateXdot, [tVec(end) tVec(1)], X0comb, options);
X_i = reshape(Xout(end,1:12), 12, 1)';
phi_i = reshape(Xout(end,13:end), 12, 12);

xhat = phi_i * xhat
P = phi_i * P * phi_i';

rmsPre = sqrt(ssrPre/length(tVec))
toc

% Plot residuals
figure(1); clf; hold on;
plot(stationDat{1}(:,1), stationDat{1}(:,2), 'b', 'linewidth', 2);
plot(stationDat{2}(:,1), stationDat{2}(:,2), 'g', 'linewidth', 2);
plot(stationDat{3}(:,1), stationDat{3}(:,2), 'r', 'linewidth', 2);
legend('Station 1', 'Station 2', 'Station 3');
xlabel('Time (s)');
ylabel('Residual (m)');

figure(2); clf; hold on;
plot(stationDatPost{1}(:,1), stationDatPost{1}(:,2), 'b', 'linewidth', 2);
plot(stationDatPost{2}(:,1), stationDatPost{2}(:,2), 'g', 'linewidth', 2);
plot(stationDatPost{3}(:,1), stationDatPost{3}(:,2), 'r', 'linewidth', 2);
legend('Station 1', 'Station 2', 'Station 3');
xlabel('Time (s)');
ylabel('Residual (m)');

    



load task_real;
IMUTs = diff(timeIMU);
dt = mean(IMUTs);
K = size(zAcc,2);
%% Measurement noise
% GNSS Position  measurement
...

% accelerometer
qA = ...^2; % accelerometer measurement noise covariance
qAb = ...^2; % accelerometer bias driving noise covariance
pAcc = ...; % accelerometer bias reciprocal time constant

qG = ...^2; % gyro measurement noise covariance
qGb = ...^2;  % gyro bias driving noise covariance
pGyro = ...; % gyrp bias reciprocal time constant


%% Estimator
eskf = ESKF(qA, qG, qAb, qGb, pAcc, pGyro);
eskf.Sa = S_a; % set the accelerometer correction matrix
eskf.Sg = S_g; % set the gyro correction matrix
%% Allocate
xest = zeros(16, K);
Pest = zeros(15, 15, K);

xpred = zeros(16, K);
Ppred = zeros(15, 15, K);

%% initialize
xpred(1:3, 1) = [0, 0, 0]'; % starting 5 meters above ground
%xpred(4:6, 1) = [0.0173; -0.0129; 0.1983]; %
xpred(7, 1) = cosd(45); % nose to east, right to south and belly down.
xpred(10, 1) = sind(45);

Ppred(1:3, 1:3, 1) = 10^2*eye(3); 
Ppred(4:6, 4:6, 1) = 3^2*eye(3);
Ppred(7:9, 7:9, 1) = (pi/30)^2 * eye(3); % error rotation vector (not quat)
Ppred(10:12, 10:12, 1) = 0.05^2 * eye(3);
Ppred(13:15, 13:15, 1) = (2e-5)^2 * eye(3);

%% run
N = K;
GNSSk = 1;
for k = 1:N
    t = timeIMU(k);
    
    if mod(k, 1000) == 0
        fprintf('time %.3f at step %d\n', t - IMUtime(1), k);
    end
    
    if timeGNSS(GNSSk) < t
        NIS(GNSSk) = ...;
        [xest(:, k), Pest(:, :, k)] = ...;
        GNSSk = GNSSk + 1;
    
        if any(any(~isfinite(Pest(:, :, k))))
            error('not finite Pest at time %d',k)
        end
    else % no updates so estimate = prediction
        xest(:, k) = ...;
        Pest(:, :, k) = ...;
    end

    if k < K
        [xpred(:, k + 1),  Ppred(:, :, k + 1)] = ...;
        
        if any(any(~isfinite(Ppred(:, :, k + 1))))
            error('not finite Ppred at time %d', k + 1)
        end
    end  
end
%% plots
figure(1);
clf;
plot3(xest(2, 1:N), xest(1, 1:N), -xest(3, 1:N));
hold on;
plot3(zGNSS(2, 1:GNSSk), zGNSS(1, 1:GNSSk), -zGNSS(3, 1:GNSSk))
grid on; axis equal
xlabel('East [m]')
ylabel('North [m]')
zlabel('Altitude [m]')

%%
eul = quat2eul(xest(7:10, :))*180/pi;
figure(2); clf; hold on;
subplot(5,1,1);
plot(timeIMU(1:N) - timeIMU(1), xest(1:3, 1:N))
grid on;
ylabel('NED position [m]')
subplot(5,1,2);
plot(timeIMU(1:N) - timeIMU(1), xest(4:6, 1:N))
grid on;
ylabel('Velocitites [m/s]')
subplot(5,1,3);
plot(timeIMU(1:N) - timeIMU(1), eul(:, 1:N))
grid on;
ylabel('euler angles [deg]')
legend('\phi', '\theta', '\psi')
subplot(5, 1, 4)
plot(timeIMU(1:N) - timeIMU(1), xest(11:13, 1:N))
grid on;
ylabel('Accl bias [m/s^2]')
subplot(5, 1, 5)
plot(timeIMU(1:N) - timeIMU(1), xest(14:16, 1:N)*180/pi * 3600)
grid on;
ylabel('Gyro bias [rad/s]')

figure(3);
alpha = 0.05;
CI3 = chi2inv([alpha/2; 1 - alpha/2; 0.5], 3);
clf;
plot(timeGNSS(1:(GNSSk - 1)) - timeIMU(1), NIS);
grid on;
hold on;
plot([0, timeIMU(N) - timeIMU(1)], (CI3*ones(1,2))', 'r--');
insideCI = mean((CI3(1) <= NIS).* (NIS <= CI3(2)));
title(sprintf('NIS (%.3g%% inside %.3g%% confidence intervall)', 100*insideCI, 100*(1 - alpha)));

figure(4); clf;
gaussCompare = sum(randn(3, numel(NIS)).^2, 1);
boxplot([NIS', gaussCompare'],'notch','on',...
        'labels',{'NIS','gauss'});
grid on;
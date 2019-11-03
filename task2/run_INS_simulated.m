load task_simulation.mat;
dt = mean(diff(timeIMU));
steps = size(zAcc,2);

%% Measurement noise
% GNSS Position  measurement
p_std = [..., ..., ...]'; % Measurement noise
RGNSS = diag(p_std.^2);

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
xest = zeros(16, steps);
Pest = zeros(15, 15, steps);

xpred = zeros(16, steps);
Ppred = zeros(15, 15, steps);

%% initialize
xpred(1:3, 1) = [0, 0, -5]'; % starting 5 meters above ground
xpred(4:6, 1) = [20, 0, 0]'; % starting at 20 m/s due north
xpred(7, 1) = 1; % no initial rotation: nose to north, right to East and belly down.

Ppred(1:3, 1:3, 1) = ...; 
Ppred(4:6, 4:6, 1) = ...;
Ppred(7:9, 7:9, 1) = ...; % error rotation vector (not quat)
Ppred(10:12, 10:12, 1) = ...;
Ppred(13:15, 13:15, 1) = ...;

%% run
N = 90000;
GNSSk = 1;
for k = 1:N
    if  timeIMU(k) >= timeGNSS(GNSSk)
        NIS(GNSSk) = ...;
        [xest(:, k), Pest(:, :, k)] = ...;
        GNSSk = GNSSk  + 1;
        
        % sanity check, remove for some minor speed
        if any(any(~isfinite(Pest(:, :, k))))
            error('not finite Pest at time %d',k)
        end
        
    else % no updates so estimate = prediction
        xest(:, k) = ... ;
        Pest(:, :, k) = ... ;
    end
    
    deltaX(:, k) = ...;
    [NEES(:, k), NEESpos(:, k), NEESvel(:, k), NEESatt(:, k), NEESaccbias(:, k), NEESgyrobias(:, k)] = ...
        ...;
    
    if k < N
        [xpred(:, k+1),  Ppred(:, :, k+1)] = ...;
        % sanity check, remove for speed
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

% state estimate plot
eul = quat2eul(xest(7:10, :));
eul_true = quat2eul(xtrue(7:10, :));
figure(2); clf; hold on;

subplot(5,1,1);
plot((0:(N-1))*dt, xest(1:3, 1:N))
grid on;
ylabel('NED position [m]')
legend('North', 'East', 'Down')

subplot(5,1,2);
plot((0:(N-1))*dt, xest(4:6, 1:N))
grid on;
ylabel('Velocitites [m/s]')
legend('North', 'East', 'Down')

subplot(5,1,3);
plot((0:(N-1))*dt, eul(:, 1:N)*180/pi)
hold on
%plot((0:(N-1))*dt, euler_out(:, 1:N)*180/pi)
grid on;
ylabel('euler angles [deg]')
legend('\phi', '\theta', '\psi')

subplot(5, 1, 4)
plot((0:(N-1))*dt, xest(11:13, 1:N))
grid on;
ylabel('Accl bias [m/s^2]')
legend('x', 'y', 'z')
subplot(5, 1, 5)
plot((0:(N-1))*dt, xest(14:16, 1:N)*180/pi * 3600)
grid on;
ylabel('Gyro bias [deg/h]')
legend('x', 'y', 'z')

suptitle('States estimates');

% state error plots
figure(3); clf; hold on;

subplot(5,1,1);
plot((0:(N-1))*dt, deltaX(1:3,:))
grid on;
ylabel('NED position error [m]')
legend(sprintf('North (%.3g)', sqrt(mean(deltaX(1, 1:N).^2))),...
    sprintf('East (%.3g)', sqrt(mean(deltaX(2, 1:N).^2))),...
    sprintf('Down (%.3g)', sqrt(mean(deltaX(3, 1:N).^2))))

subplot(5,1,2);
plot((0:(N-1))*dt, deltaX(4:6, 1:N))
grid on;
ylabel('Velocitites error [m/s]')
legend(sprintf('North (%.3g)', sqrt(mean(deltaX(4, 1:N).^2))),...
    sprintf('East (%.3g)', sqrt(mean(deltaX(5, 1:N).^2))),...
    sprintf('Down (%.3g)', sqrt(mean(deltaX(6, 1:N).^2))))

subplot(5,1,3);
plot((0:(N-1))*dt, wrapToPi(eul(:, 1:N) - eul_true(:, 1:N))*180/pi)
grid on;
ylabel('euler angles error [deg]')
legend(sprintf('\\phi (%.3g)', sqrt(mean((eul(1, 1:N) - eul_true(1, 1:N)).^2))),...
    sprintf('\\theta (%.3g)', sqrt(mean((eul(2, 1:N) - eul_true(2, 1:N)).^2))),...
    sprintf('\\psi (%.3g)', sqrt(mean((eul(3, 1:N) - eul_true(3, 1:N)).^2))))

subplot(5, 1, 4)
plot((0:(N-1))*dt, deltaX(10:12, 1:N))
grid on;
ylabel('Accl bias error [m/s^2]')
legend(sprintf('x (%.3g)', sqrt(mean(deltaX(1, 1:N).^2))),...
    sprintf('y (%.3g)', sqrt(mean(deltaX(11, 1:N).^2))),...
    sprintf('z (%.3g)', sqrt(mean(deltaX(12, 1:N).^2))))

subplot(5, 1, 5)
plot((0:(N-1))*dt, deltaX(13:15, 1:N)*180/pi);
grid on;
ylabel('Gyro bias error [deg/s]')
legend(sprintf('x (%.3g)', sqrt(mean(((deltaX(13, 1:N))*180/pi).^2))),...
    sprintf('y (%.3g)', sqrt(mean(((deltaX(14, 1:N))*180/pi).^2))),...
    sprintf('z (%.3g)', sqrt(mean(((deltaX(15, 1:N))*180/pi).^2))))

suptitle('States estimate errors');

% error distance plot
figure(4); clf; hold on;
subplot(2,1,1); hold on;
plot((0:(N-1))*dt, sqrt(sum(deltaX(1:3, 1:N).^2,1)))
plot((0:100:(N-1))*dt, sqrt(sum((xtrue(1:3, 100:100:N) - zGNSS(:, 1:GNSSk)).^2,1)))
ylabel('Position error [m]')
grid on;
legend(sprintf('estimation error (%.3g)',sqrt(mean(sum(deltaX(1:3, 1:N).^2,1))) ),...
    sprintf('measurement error (%.3g)', sqrt(mean(sum((xtrue(1:3, 100:100:N) - zGNSS(:, 1:GNSSk)).^2,1)))));

subplot(2,1,2);
plot((0:(N-1))*dt, sqrt(sum(deltaX(4:6, 1:N).^2, 1)))
ylabel('Speed error [m/s]');
title(sprintf('RMSE: %.3g', sqrt(mean(sum(deltaX(4:6, 1:N).^2, 1)))));
grid on;

%% CONSISTENCY
alpha = 0.05;
CI15 = chi2inv([alpha/2; 1 - alpha/2; 0.5], 15);
CI3 = chi2inv([alpha/2; 1 - alpha/2; 0.5], 3);

figure(5); clf;
subplot(7,1,1);
plot((0:(N-1))*dt, NEES);
grid on;
hold on;
plot([0, N-1]*dt, (CI15*ones(1,2))', 'r--');
insideCI = mean((CI15(1) <= NEES).* (NEES <= CI15(2)));
title(sprintf('total NEES (%.3g%% inside %.3g%% confidence intervall)', 100*insideCI, 100*(1 - alpha)));

subplot(7,1,2);
plot((0:(N-1))*dt, NEESpos);
grid on;
hold on;
plot([0, N-1]*dt, (CI3*ones(1,2))', 'r--');
insideCI = mean((CI3(1) <= NEESpos).* (NEESpos <= CI3(2)));
title(sprintf('position NEES (%.3g%% inside %.3g%% confidence intervall)', 100*insideCI, 100*(1 - alpha)));

subplot(7,1,3);
plot((0:(N-1))*dt, NEESvel);
grid on;
hold on;
plot([0, N-1]*dt, (CI3*ones(1,2))', 'r--');
insideCI = mean((CI3(1) <= NEESvel).* (NEESvel <= CI3(2)));
title(sprintf('velocity NEES (%.3g%% inside %.3g%% confidence intervall)', 100*insideCI, 100*(1 - alpha)));

subplot(7,1,4);
plot((0:(N-1))*dt, NEESatt);
grid on;
hold on;
plot([0, N-1]*dt, (CI3*ones(1,2))', 'r--');
insideCI = mean((CI3(1) <= NEESatt).* (NEESatt <= CI3(2)));
title(sprintf('attitude NEES (%.3g%% inside %.3g%% confidence intervall)', 100*insideCI, 100*(1 - alpha)));

subplot(7,1,5);
plot((0:(N-1))*dt, NEESaccbias);
grid on;
hold on;
plot([0, N-1]*dt, (CI3*ones(1,2))', 'r--');
insideCI = mean((CI3(1) <= NEESaccbias).* (NEESaccbias <= CI3(2)));
title(sprintf('accelerometer bias NEES (%.3g%% inside %.3g%% confidence intervall)', 100*insideCI, 100*(1 - alpha)));

subplot(7,1,6);
plot((0:(N-1))*dt, NEESgyrobias);
grid on;
hold on;
plot([0, N-1]*dt, (CI3*ones(1,2))', 'r--');
insideCI = mean((CI3(1) <= NEESgyrobias).* (NEESgyrobias <= CI3(2)));
title(sprintf('gyro bias NEES (%.3g%% inside %.3g%% confidence intervall)', 100*insideCI, 100*(1 - alpha)));

subplot(7,1,7)
plot(0:(numel(NIS)-1), NIS);
grid on;
hold on;
plot([0, N-1]*dt, (CI3*ones(1,2))', 'r--');
insideCI = mean((CI3(1) <= NIS).* (NIS <= CI3(2)));
title(sprintf('NIS (%.3g%% inside %.3g%% confidence intervall)', 100*insideCI, 100*(1 - alpha)));

% boxplot
figure(6)
subplot(1,3,1)
gaussCompare = sum(randn(3, numel(NIS)).^2, 1);
boxplot([NIS', gaussCompare'],'notch','on',...
        'labels',{'NIS','gauss'})
grid on
subplot(1,3,2)
gaussCompare15 = sum(randn(15, N).^2, 1);
gaussCompare3 = sum(randn(3, N).^2, 1);
boxplot([NEES', gaussCompare15'],'notch', 'on', 'labels',{'NEES','gauss(15dim)'});
grid on;
subplot(1,3,3)
boxplot([NEESpos', NEESvel', NEESatt', NEESaccbias', NEESgyrobias', gaussCompare3'],...
    'notch', 'on', 'labels',{'NEESpos', 'NEESvel', 'NEESatt', 'NEESaccbias', 'NEESgyrobias', 'gauss(3dim)'})
grid on;
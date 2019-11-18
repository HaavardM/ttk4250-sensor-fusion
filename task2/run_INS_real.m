%% Parameters
showplt_track = true;
showplt_estimates = true;
showplt_nis = true;
showplt_boxplot = true;
showplt_nis_colored_track = true; % Must be true to print
showplt_ppred_colored_track = true; % Must be true to print
showplt_covariance = true;

%% Load data
load task_real;
IMUTs = diff(timeIMU);
dt = mean(IMUTs);
K = size(zAcc,2);
%% Measurement noise

% accelerometer
p_std = (mean(GNSSaccuracy, 2)) * [1, 1 , 1]'; % Measurement noise
%RGNSS = diag(p_std.^2);
RGNSS = @(k) diag(((0.20*GNSSaccuracy(k))^2)*[1 1 1]);

% accelerometer
qA = (1.167e-3)^2;% accelerometer measurement noise covariance
qAb = (1.5e-3)^2; % accelerometer bias driving noise covariance
pAcc = 1e-8; % accelerometer bias reciprocal time constant

qG = (deg2rad(2.5e-3))^2; % gyro measurement noise covariance
qGb = (8e-6)^2;  % gyro bias driving noise covariance
pGyro = 1e-8; % gyrp bias reciprocal time constant



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

Ppred(1:3, 1:3, 1) = (2e-1)^2*eye(3); % pos
Ppred(4:6, 4:6, 1) = (3e-4)^2*eye(3); % vel
Ppred(7:9, 7:9, 1) = (1e-1*(pi/30))^2 * eye(3); % error rotation (vector (not quat)
Ppred(10:12, 10:12, 1) = 0.02^2 * eye(3); % acc bias
Ppred(13:15, 13:15, 1) = (1e-4)^2 * eye(3); % gyro bias

NIS = zeros(1, 3682);

%% run
N = K;
GNSSk = 1;
for k = 1:N
    t = timeIMU(k);
    if mod(k, 1000) == 0
        fprintf('time %.3f at step %d\n', t - IMUTs(1), k);
    end
    
    if timeGNSS(GNSSk) < t
        NIS(GNSSk) = eskf.NISGNSS(xpred(:, k), Ppred(:, :, k), zGNSS(:, GNSSk), RGNSS(GNSSk), leverarm);
        [xest(:, k), Pest(:, :, k)] = eskf.updateGNSS(xpred(:, k), Ppred(:, :, k), zGNSS(:, GNSSk), RGNSS(GNSSk));
        GNSSk = GNSSk + 1;
    
        if any(any(~isfinite(Pest(:, :, k))))
            error('not finite Pest at time %d',k)
        end
    else % no updates so estimate = prediction
        xest(:, k) = xpred(:, k) ;
        Pest(:, :, k) = Ppred(:, :, k) ;
    end

    if k < K
        [xpred(:, k+1),  Ppred(:, :, k+1)] = eskf.predict(xest(:, k), Pest(:, :, k), zAcc(:, k+1), zGyro(:, k+1), IMUTs(k));
        
        % Sanity check: Remove for speeeeeeeeeeeeeeeeed
%         if any(any(~isfinite(Ppred(:, :, k + 1))))
%             error('not finite Ppred at time %d', k + 1)
%         end
    end  
end
%% plots
if exist('showplt_track') && showplt_track
    figure(1);
else
    figure("visible", "off");
end
clf;
plot3(xest(2, 1:N), xest(1, 1:N), -xest(3, 1:N));
hold on;
plot3(zGNSS(2, 1:GNSSk), zGNSS(1, 1:GNSSk), -zGNSS(3, 1:GNSSk))
grid on; axis equal
xlabel('East [m]')
ylabel('North [m]')
zlabel('Altitude [m]')
printplot(gcf, "a2-real-track.pdf");

%%
eul = zeros(3, N);
for i =1:N
   eul(:, i) = quat2eul(xest(7:10, i));
end
%% Plot stuff
if exist('showplt_estimates') && showplt_estimates
    figure(2); clf; hold on;
else
    figure('visible', 'off'); clf; hold on;
end
subplot(5,1,1);
plot(timeIMU(1:N) - timeIMU(1), xest(1:3, 1:N)); hold on;
plot(timeGNSS(1:GNSSk) - timeGNSS(1), zGNSS(1:3, 1:GNSSk));
grid on;
ylabel('NED position [m]')
legend();
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
printplot(gcf, "a2-real-estimates.pdf");

if exist('showplt_nis') && showplt_nis
    fig3 = figure(3);
else
    fig3 = figure('visible', 'off');
end
alpha = 0.05;
CI3 = chi2inv([alpha/2; 1 - alpha/2; 0.5], 3);
clf;
plot(timeGNSS(1:(GNSSk - 1)) - timeIMU(1), NIS);
grid on;
hold on;
plot([0, timeIMU(N) - timeIMU(1)], (CI3*ones(1,2))', 'r--');
insideCI = mean((CI3(1) <= NIS).* (NIS <= CI3(2)));
title(sprintf('NIS (%.3g%% inside %.3g%% confidence intervall)', 100*insideCI, 100*(1 - alpha)));
printplot(fig3, "a2-real-nis.pdf");

if exist('showplt_boxplot') && showplt_boxplot
    figure(4); clf;
else
    figure('visible', 'off'); clf;
end
gaussCompare = sum(randn(3, numel(NIS)).^2, 1);
boxplot([NIS', gaussCompare'],'notch','on',...
        'labels',{'NIS','gauss'});
grid on;
printplot(gcf, "a2-real-boxplot.pdf");

%% NIS colored track plot
start = 1;
if exist('showplt_nis_colored_track') && showplt_nis_colored_track
    plotcoloredtrack(zGNSS(:, start:GNSSk-1), NIS(start:GNSSk-1), "NIS colored track (xy projection)", 5, 5);
    printplot(gcf, "a2-real-nis-colored-track.pdf");
end

%%
states(1) = "position";
states(4) = "velocity";
states(7) = "attitude";
states(10) = "acc bias";
states(13) = "gyro bias";

%% Ppred colored track
if exist('showplt_ppred_colored_track') && showplt_ppred_colored_track
    plot_no = 0;
    for state = 1:3:13
        plot_no = plot_no + 1;
        k = 1;
        for i = 1:N
            if timeIMU(i) > timeGNSS(k)
                GNSS_Ppred_norm(:, :, k) = norm(Ppred(state:state+2, state:state+2, i));
                k = k + 1;
            end
        end
        start = 1;
        plotcoloredtrack(zGNSS(:, start:GNSSk-1), GNSS_Ppred_norm(:,start:GNSSk-1), sprintf("%s Ppred colored track", states(state)), 3000, 6 + state);
        printplot(gcf, sprintf("a2-real-%s-ppred-colored-track.pdf", states(state)));
    end
end
%% Covariance plot
if exist('showplt_covariance') && showplt_covariance
    fig9 = figure(90);
else
    fig9 = figure('visible', 'off');
end
for state = 1:3:13
    subplot(5,1,plot_no);
    plot(GNSS_Ppred_norm(:, start:GNSSk-1));
    title(sprintf("%s Ppred", states(state)));
end
printplot(gcf, "a2-real-covariance.pdf");
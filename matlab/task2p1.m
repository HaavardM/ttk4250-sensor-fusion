% load data
usePregen = false; % you can generate your own data if set to false
if usePregen
    load task2data.mat;
else
    K = 100;
    Ts = 2.5;
    r = 5;
    q = [0.005, 1e-6*pi]; % q(1): CV noise (effective all the time), q(2): CT noise (effective only in turns)
    init.x = [0; 0; 2; 0; 0];
    init.P = diag([25, 25, 3, 3, 0.0005].^2);
    % detection and false alarm
    PDtrue = 0.9;
    lambdatrue = 1e-4;
    [Xgt, Z, a] = simulate_atc_track(Ts, K, q, r, init, PDtrue, lambdatrue, false);
end


% plot measurements close to the trajectory
figure(1); clf; hold on; grid on;
Zplotdata = [];
plotMeasDist = 45;
for k = 1:K
   toPlot = false(size(Z{k},2),1);
   for j = 1:size(Z{k}, 2)
        v = Z{k}(:, j) - Xgt(1:2, k);
        toPlot(j) = v' * v <= plotMeasDist^2;
   end
   Zplotdata = [Zplotdata, Z{k}(:, toPlot)];
end
set(gca, 'ColorOrderIndex', 2)
scatter(Zplotdata(1,:), Zplotdata(2,:));
set(gca, 'ColorOrderIndex', 1)
plot(Xgt(1,:),Xgt(2,:), 'LineWidth',1.5);
title('True trajectory and the nearby measurements')

% play measurement movie. Remember that you can click the stop button.
% figure(2); clf; grid on;
% set(gcf,'Visible','on')
% plotpause = 0.1; % sets a pause in between time steps if it goes to fast
% minZ = min(Xgt(1:2, :), [], 2);
% maxZ = max(Xgt(1:2, :), [], 2);
% minZ = minZ - 0.05*(maxZ - minZ);
% maxZ = maxZ + 0.05*(maxZ - minZ);
% for k = 1:K
%     scatter(Z{k}(1,:), Z{k}(2,:))
%     title(sprintf('measurements at step %d', k))
%     axis([minZ(1), maxZ(1), minZ(2), maxZ(2)])
%     grid on;
%     drawnow;
%     pause(plotpause);
% end

% Look at individual EKF-PDAs
r = 10;
qCV = 1.0;
qCT = [6.5, 6.5];

lambda = 1e-3;
PD = 0.8;
gateSize = 5^2;
% choose model to tune
s = 1;

% make models
models =  cell(2,1);
models{1} = EKF(discreteCVmodel(qCV, r));
models{2} = EKF(discreteCTmodel(qCT, r));

tracker =  cell(2,1);
tracker{1} = PDAF(models{1}, lambda, PD, gateSize);
tracker{2} = PDAF(models{2}, lambda, PD, gateSize);

% % % % allocate
xbar = zeros(5, K);
Pbar = zeros(5, 5, K);
xhat = zeros(5, K);
Phat = zeros(5, 5, K);
NEES = zeros(K, 1);
NEESpos = zeros(K, 1);
NEESvel = zeros(K, 1);

% initialize filter
xbar(:, 1) = [0; 0; 2; 0; 0];
Pbar(:, : ,1) = diag([25, 25, 3, 3, 0.0005].^2);

% filter
for k = 1:K
    [xhat(: , k) , Phat(: , :, k)] = tracker{s}.update(Z{k} , xbar(: , k), Pbar(: , :, k));
    NEES(k) = ((xhat(1:4, k) - Xgt(1:4, k))'/ squeeze(Phat(1:4, 1:4, k))) * (xhat(1:4, k) - Xgt(1:4, k));
    NEESpos(k) = ((xhat(1:2, k) - Xgt(1:2, k))'/ squeeze(Phat(1:2, 1:2, k))) * (xhat(1:2, k) - Xgt(1:2, k));
    NEESvel(k) = ((xhat(3:4, k) - Xgt(3:4, k))'/ squeeze(Phat(3:4, 3:4, k))) * (xhat(3:4, k) - Xgt(3:4, k));
    if k < K
        [xbar(:, k+1), Pbar(:, :, k+1)] = tracker{s}.predict(xhat(:, k), Phat(:, :, k), Ts);
    end
end

% errors
poserr = sqrt(sum((xhat(1:2,:) - Xgt(1:2,:)).^2, 1));
posRMSE = sqrt(mean(poserr.^2));
velerr = sqrt(sum((xhat(3:4, :) - Xgt(3:4, :)).^2, 1));
velRMSE = sqrt(mean(velerr.^2));

% consistency
chi2inv([0.025, 0.975], K*2)/K
ANEESpos = mean(NEESpos)
ANEESvel = mean(NEESvel)

chi2inv([0.025, 0.975], K*4)/K
ANEES = mean(NEES)

% plot
figure(3); clf; hold on; grid on;
plot(xhat(1,:), xhat(2,:));
plot(Xgt(1, :), Xgt(2, :))
axis('equal')
%scatter(Z(1,:), Z(2, :));
title(sprintf('posRMSE = %.3f, velRMSE = %.3f',posRMSE, velRMSE))

figure(4); clf; hold on; grid on;
plot(xhat(5,:))
plot(Xgt(5,:))
ylabel('omega')

figure(5); clf;
subplot(3,1,1);
plot(NEES); grid on; hold on;
ylabel('NEES');
ciNEES = chi2inv([0.05, 0.95], 4);
inCI = sum((NEES >= ciNEES(1)) .* (NEES <= ciNEES(2)))/K * 100;
plot([1,K], repmat(ciNEES',[1,2])','r--')
text(104, -5, sprintf('%.2f%% inside CI', inCI),'Rotation',90);

subplot(3,1,2);
plot(NEESpos); grid on; hold on;
ylabel('NEESpos');
ciNEES = chi2inv([0.05, 0.95], 2);
inCI = sum((NEESpos >= ciNEES(1)) .* (NEESpos <= ciNEES(2)))/K * 100;
plot([1,K], repmat(ciNEES',[1,2])','r--')
text(104, -5, sprintf('%.2f%% inside CI', inCI),'Rotation',90);

subplot(3,1,3);
plot(NEESvel); grid on; hold on;
ylabel('NEESvel');
ciNEES = chi2inv([0.05, 0.95], 2);
inCI = sum((NEESvel >= ciNEES(1)) .* (NEESvel <= ciNEES(2)))/K * 100;
plot([1,K], repmat(ciNEES',[1,2])','r--')
text(104, -5, sprintf('%.2f%% inside CI', inCI),'Rotation',90);

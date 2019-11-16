% load data
load joyridedata.mat;

addpath('./plotfuncs');
%% Parameters
fignum = 1;

r = 200;
lambda = 5e-5;
PD = 0.85;
gateSize = 5^2;

% dynamic models
qCV = 4;
qCT = [5, 0.00002];
qCVh = 100;
modIdx = 1:3; 
M = numel(modIdx);
PI = [0.95    0.04    0.25
      0.04    0.95    0.25
      0.01    0.01    0.50];
  
PI = PI(modIdx, modIdx); % select the models to use
PI = PI./sum(PI,1); % be sure to normalize
assert(all(sum(PI, 1) - 1 < eps),'columns of PI must sum to 1');

x0 = [7100; 3630; 0; 0; 0]; % taken from gt
P0 = diag([25, 25, 10, 10, pi/6].^2); % seems reasonable?

%% plot measurements close to the trajectory
figure(fignum); clf; hold on; grid on;
fignum = fignum + 1;
Zmat = cell2mat(Z');
Zplotdata = [];
plotMeasDist = 200;
for k = 1:K
   toPlot = false(size(Z{k},2),1);
   closecount(k) = 0;
   for j = 1:size(Z{k}, 2)
        v = Z{k}(:, j) - Xgt(1:2, k);
        toPlot(j) = v' * v <= plotMeasDist^2;
        if toPlot(j) 
            closecount(k) = closecount(k) + 1;
        end
   end
   Zplotdata = [Zplotdata, Z{k}(:, toPlot)];
end
set(gca, 'ColorOrderIndex', 2)
scatter(Zplotdata(1,:), Zplotdata(2,:));
set(gca, 'ColorOrderIndex', 1)
plot(Xgt(1,:),Xgt(2,:), 'LineWidth',1.5);
title('True trajectory and the nearby measurements')

% some minor data analysis for the interested person
i = 0;
maxdist = 30;
for k = 1:K
    if size(Z{k}, 2) > 0
        vall = Xgt(1:2,k) - Z{k};
        dsquared = sum((vall).^2, 1);
        closecount(k) = sum(dsquared <= maxdist^2);
        [mindist, ind] = min(dsquared);
        if mindist < maxdist^2
            i = i + 1;
            d(i) = mindist;
            v(:, i) = vall(:, ind);
        end
    end
end
stepsWithClose = i;
meanMeasError = mean(v, 2);
covMeasError = cov(v');
stds = sqrt(diag(cov(v')));
covEllSize = sqrt(eig(cov(v')));
meanMinMeasDists = mean(d);
meanNumberOfCloseMeasurements = mean(closecount);

% make model
tracker =  cell(2,1);
tracker{1} = PDAF(EKF(discreteCVmodel(qCV, r)), lambda, PD, gateSize);
tracker{2} = PDAF(EKF(discreteCTmodel(qCT, r)), lambda, PD, gateSize);

% allocate
xbar = zeros(5, K);
Pbar = zeros(5, 5, K);
xest = zeros(5, K);
Phat = zeros(5, 5, K);
NEES = zeros(K, 1);
NEESpos = zeros(K, 1);
NEESvel = zeros(K, 1);

s = 2;
% initialize
xbar(:, 1) = x0; % simply same for all modes
Pbar(:, :, 1) = P0; %simply same for all modes

% filter
for k=1:K
    [xest(:, k), Phat(:, :, k)] = tracker{s}.update(Z{k}, xbar(:, k), Pbar(:, :, k));
    NEES(k) =  (xest(1:4, k) - Xgt(1:4, k))' * (Phat(1:4, 1:4, k )\ (xest(1:4, k) - Xgt(1:4, k)));
    NEESpos(k) = (xest(1:2, k) - Xgt(1:2, k))' * (Phat(1:2, 1:2, k )\ (xest(1:2, k) - Xgt(1:2, k)));
    NEESvel(k) = (xest(3:4, k) - Xgt(3:4, k))' * (Phat(3:4, 3:4, k )\ (xest(3:4, k) - Xgt(3:4, k)));
    if k < K
        [xbar(:, k+1), Pbar(:, :,k+1)] = tracker{s}.predict(xest(:, k), Phat(:, :,k), Ts(k));
    end
end

% errors
poserr = sqrt(sum((xest(1:2,:) - Xgt(1:2,:)).^2, 1));
posRMSE = sqrt(mean(poserr.^2)); % not true RMSE (which is over monte carlo simulations)
velerr = sqrt(sum((xest(3:4, :) - Xgt(3:4, :)).^2, 1));
velRMSE = sqrt(mean(velerr.^2)); % not true RMSE (which is over monte carlo simulations)
peakPosDeviation = max(poserr);
peakVelDeviation = max(velerr);

% consistency
CI2K = chi2inv([0.025, 0.975], K*2)/K;
ANEESpos = mean(NEESpos);
ANEESvel = mean(NEESvel);

CI4K = chi2inv([0.025, 0.975], K*4)/K;
ANEES = mean(NEES);

% plot
figure(fignum); clf; hold on; grid on;
fignum = fignum + 1;
plot(xest(1,:), xest(2,:));
plot(Xgt(1,:), Xgt(2, :));
scatter(Zmat(1,:), Zmat(2,:))
axis('equal')
title(sprintf('posRMSE = %.3f, velRMSE = %.3f, peakPosDev = %.3f, peakVelDev = %.3f',posRMSE, velRMSE, peakPosDeviation, peakVelDeviation))
print -depsc plots/task3/a1_task3_tracked_path

figure(fignum); clf; hold on; grid on;
fignum = fignum + 1;
subplot(3,1,1);
hold on; grid on;
plot(sqrt(sum(xest(3:4,:).^2,1)))
plot(sqrt(sum(Xgt(3:4,:).^2,1)))
ylabel('speed')
subplot(3,1,2);
hold on; grid on;
plot(atan2(xest(4,:), xest(3,:)))
plot(atan2(Xgt(4,:), Xgt(3,:)))
ylabel('theta')
subplot(3,1,3)
hold on; grid on;
plot(diff(unwrap(atan2(xest(4,:), xest(3,:))))./Ts')
plot(diff(unwrap(atan2(Xgt(4,:), Xgt(3,:))))./Ts')
ylabel('omega')
print -depsc plots/task3/a1_task3_states

figure(fignum); clf;
fignum = fignum + 1;
subplot(2,1,1); 
plot(poserr); grid on;
ylabel('position error')
subplot(2,1,2);
plot(velerr); grid on;
ylabel('velocity error')
print -depsc plots/task3/a1_task3_error

figure(fignum); clf;
fignum = fignum + 1;
subplot(3,1,1);
plot(NEES); grid on; hold on;
ylabel('NEES');
ciNEES = chi2inv([0.05, 0.95], 4);
inCI = sum((NEES >= ciNEES(1)) .* (NEES <= ciNEES(2)))/K * 100;
plot([1,K], repmat(ciNEES',[1,2])','r--')
text(K*1.04, -5, sprintf('%.2f%% inside CI', inCI),'Rotation',90);


subplot(3,1,2);
plot(NEESpos); grid on; hold on;
ylabel('NEESpos');
ciNEES = chi2inv([0.05, 0.95], 2);
inCI = sum((NEESpos >= ciNEES(1)) .* (NEESpos <= ciNEES(2)))/K * 100;
plot([1,K], repmat(ciNEES',[1,2])','r--')
text(K*1.04, -5, sprintf('%.2f%% inside CI', inCI),'Rotation',90);

subplot(3,1,3);
plot(NEESvel); grid on; hold on;
ylabel('NEESvel');
ciNEES = chi2inv([0.05, 0.95], 2);
inCI = sum((NEESvel >= ciNEES(1)) .* (NEESvel <= ciNEES(2)))/K * 100;
plot([1,K], repmat(ciNEES',[1,2])','r--')
text(K*1.04, -5, sprintf('%.2f%% inside CI', inCI),'Rotation',90);
print -depsc plots/task3/a1_task3_NEES


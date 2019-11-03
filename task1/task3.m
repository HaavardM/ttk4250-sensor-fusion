% load data
load joyridedata.mat;

addpath('./plotfuncs');
%% Parameters
fignum = 1;

r = 6;
lambda = 1e-4;
PD = 0.95;
gateSize = 5^2;

% dynamic models
qCV = 3;
qCT = [5, 0.000005];
qCVh = 50;
modIdx = 1:3; 
M = numel(modIdx);
PI = [0.95    0.025    0.025
      0.025    0.95    0.025
      0.025    0.025    0.95];
  
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
xbar(:, 1) = [7100; 3630; -6; -3; 0]; % taken from ground truth
Pbar(:, : ,1) = diag(4*[10; 10; 1; 1; pi/10].^2); % seems reasonable?

% errors
poserr = sqrt(sum((xhat(1:2,:) - Xgt(1:2,:)).^2, 1));
posRMSE = sqrt(mean(poserr.^2));
velerr = sqrt(sum((xhat(3:4, :) - Xgt(3:4, :)).^2, 1));
velRMSE = sqrt(mean(velerr.^2));


sprobs0 = [0.5; 0.4; 0.1]; assert(sum(sprobs0) == 1, 'initial mode probabilities must sum to 1');
sprobs0 = sprobs0(modIdx)/sum(sprobs0(modIdx)); % select models and normalize
assert(all(sprobs0 > 0), 'probabilities must be positive');


% make model
models =  cell(3,1);
models{1} = EKF(discreteCVmodel(qCV, r));
models{2} = EKF(discreteCTmodel(qCT, r));
models{3} = EKF(discreteCVmodel(qCVh, r));
imm = IMM(models(modIdx), PI);
tracker = IMMPDAF(imm, lambda, PD, gateSize);

% allocate
xbar = zeros(5, M, K);
Pbar = zeros(5, 5, M, K);
probbar = zeros(M, K);
xhat = zeros(5, M, K);
xest = zeros(5, K);
Pest = zeros(5, 5, K);
Phat = zeros(5, 5, M, K);
probhat = zeros(M, K);
NEES = zeros(K, 1);
NEESpos = zeros(K, 1);
NEESvel = zeros(K, 1);

% initialize
xbar(:, :, 1) = repmat(x0, [1, M]); % simply same for all modes
Pbar(:, : ,:, 1) = repmat(P0,[1,1,M]); %simply same for all modes
probbar(:, 1) = sprobs0;

% filter
for k=1:K
    [probhat(:, k), xhat(:, :, k), Phat(:, :, :, k)] = tracker.update(Z{k}, probbar(:, k), xbar(:, :, k), Pbar(:, :, :, k));
    [xest(:, k), Pest(:, :, k)] =  tracker.imm.estimate(probhat(: , k), xhat(:, :, k ) , Phat(:, :, :, k)); %... total state mean and cov
    NEES(k) =  (xest(1:4, k) - Xgt(1:4, k))' * (Pest(1:4, 1:4, k )\ (xest(1:4, k) - Xgt(1:4, k)));
    NEESpos(k) = (xest(1:2, k) - Xgt(1:2, k))' * (Pest(1:2, 1:2, k )\ (xest(1:2, k) - Xgt(1:2, k)));
    NEESvel(k) = (xest(3:4, k) - Xgt(3:4, k))' * (Pest(3:4, 3:4, k )\ (xest(3:4, k) - Xgt(3:4, k)));
    if k < K
        [probbar(:, k+1), xbar(:, :, k+1), Pbar(:, :, :, k+1)] = tracker.predict(probhat(:, k), xhat(:, :, k), Phat(:, :, :, k), Ts(k));
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

figure(fignum); clf;
fignum = fignum + 1;
plot(probhat');
legend("CV", "CT", "CVHigh");
grid on;

figure(fignum); clf;
fignum = fignum + 1;
subplot(2,1,1); 
plot(poserr); grid on;
ylabel('position error')
subplot(2,1,2);
plot(velerr); grid on;
ylabel('velocity error')

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

plotallmodecoloredtracks(xest, probhat, fignum);
fignum = fignum + 1;

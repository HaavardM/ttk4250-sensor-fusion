% Look at individual EKF-PDAs
r = 5;
qCV = 0.0079;
qCT = [0.02, 0.0005];

lambda = 1e-4;
PD = 0.95;
gateSize = 5^2;
% choose model to tune
s = 2;

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

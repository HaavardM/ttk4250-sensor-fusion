clear all; close all;
%% Load/generate data
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
    [Xgt, Z, a] = simulate_atc_track(Ts, K, q, r, init, PDtrue, lambdatrue, false, true, true);
end


%% plot measurements close to the trajectory
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

%% Task 2p1 (PDAF)
%task2p1;


%% Task2p2 (IMM-PDAF)
task2p2;

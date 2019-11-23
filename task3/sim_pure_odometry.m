%% Clear
clear all; close all;

%% Load
load simulatedSLAM;
K = numel(z);

%% Meta parameters
results_plot_mask = [true, true];

%% Allocate
xpred = cell(1, K); % Cell array of etas
xpred{1} = poseGT(:, 1);

%% Run
slam = EKFSLAM(eye(3), eye(2));
for k = 1:K
    [xpred{k+1}, ~] = slam.predict(xpred{k}, eye(3), odometry(:, k));
end

%% Plot
if results_plot_mask(1)
    fig = figure(3);
else
    fig = figure("visible", "off");
end
if any(results_plot_mask)
    k = K;
    clf;
    hold on;

    scatter(landmarks(1,:), landmarks(2,:), 'r^')

    lh1 = plot(poseGT(1, 1:k), poseGT(2,1:k), 'r', 'DisplayName', 'gt');
    lh2 = plot(cellfun(@(x) x(1), xpred), cellfun(@(x) x(2), xpred), 'b', 'DisplayName', 'est');

    axis equal;
    title('Pure odometry results')
    legend([lh1, lh2])
    grid on;
    if results_plot_mask(2)
        printplot(fig, "a3-sim-pure-odometry.pdf");
    end
end

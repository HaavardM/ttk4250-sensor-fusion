% IMM-PDA
% sensor
r = 5;
lambda = 1e-3;
PD = 0.8;
gateSize = 5^2;

% dynamic models
qCV = 0.0025;
qCT = [0.005, 0.0005];
x0 = [0; 0; 2; 0; 0];
P0 = diag([25, 25, 3, 3, 0.0005].^2);

% markov chain (you are free to parametrize this in another way)
PI11 = 0.95;
PI22 = 0.95;
p10 = 0.5;

PI = [PI11, (1 - PI22); (1 - PI11), PI22]; assert(all(sum(PI, 1) == [1, 1]),'columns of PI must sum to 1')
sprobs0 = [p10; (1 - p10)]; assert(sum(sprobs0) == 1, 'initial mode probabilities must sum to 1');


% make model
models =  cell(2,1);
models{1} = EKF(discreteCVmodel(qCV, r));
models{2} = EKF(discreteCTmodel(qCT, r));
imm = IMM(models, PI);
tracker = IMMPDAF(imm, lambda, PD, gateSize);

% allocate
xbar = zeros(5, 2, K);
Pbar = zeros(5, 5, 2, K);
probbar = zeros(2, K);
xhat = zeros(5, 2, K);
xest = zeros(5, K);
Pest = zeros(5, 5, K);
Phat = zeros(5, 5, 2, K);
probhat = zeros(2, K);
NEES = zeros(K, 1);
NEESpos = zeros(K, 1);
NEESvel = zeros(K, 1);

% initialize
xbar(:, :, 1) = repmat(x0, [1, 2]);
Pbar(:, : ,:, 1) = repmat(P0,[1,1,2]);
probbar(:, 1) = sprobs0;

% filter
for k=1:K
    [probhat(:, k), xhat(:, :, k), Phat(:, :, :, k)] = tracker.update(Z{k}, probbar(:, k), xbar(:, :, k), Pbar(:, :, :, k));
    [xest(:, k), Pest(:, :, k)] =  tracker.imm.estimate(probhat(: , k), xhat(:, :, k ) , Phat(:, :, :, k)); %... total state mean and cov
    NEES(k) =  (xest(1:4, k) - Xgt(1:4, k))' * (Pest(1:4, 1:4, k )\ (xest(1:4, k) - Xgt(1:4, k)));
    NEESpos(k) = (xest(1:2, k) - Xgt(1:2, k))' * (Pest(1:2, 1:2, k )\ (xest(1:2, k) - Xgt(1:2, k)));
    NEESvel(k) = (xest(3:4, k) - Xgt(3:4, k))' * (Pest(3:4, 3:4, k )\ (xest(3:4, k) - Xgt(3:4, k)));
    if k < K
        [probbar(:, k+1), xbar(:, :, k+1), Pbar(:, :, :, k+1)] = tracker.predict(probhat(:, k), xhat(:, :, k), Phat(:, :, :, k), Ts);
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
if plot_tracked_path
    figure(6); clf;
    subplot(2, 1, 1);
    plot(xest(1,:), xest(2,:)); hold on; grid on;
    plot(Xgt(1,:), Xgt(2, :));
    axis('equal')
    title(sprintf('posRMSE = %.3f, velRMSE = %.3f, peakPosDev = %.3f, peakVelDev = %.3f',posRMSE, velRMSE, peakPosDeviation, peakVelDeviation))
    subplot(2, 1, 2);
    plot(xest(5,:)); hold on; grid on;
    plot(Xgt(5,:));
end

if plot_mode_prob_path
    plotmodecoloredtrack(xest, probhat, 1);
end



if plot_mode_prob
    figure(8); clf;
    plot(probhat');
    grid on;
end

if plot_error
    figure(9); clf;
    subplot(2,1,1);
    plot(poserr); grid on;
    ylabel('position error')
    subplot(2,1,2);
    plot(velerr); grid on;
    ylabel('velocity error')
end

if plot_NEES
    figure(10); clf;
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
end

if plot_movie
    %estimation "movie"
    mTL = 0.2; % maximum transparancy (between 0 and 1);
    plotpause = 1; % lenght to pause between time steps;
    plotRange = 2:2; %2:K % the range to go through
    N = 50; % number of points to use for ellipse;

    %k = 31; assert(all([k > 1, k <= K]), 'K must be in proper range')
    figure(11); clf; grid on; hold on; axis equal;
    set(gcf,'Visible','on')
    co = get(gca, 'colorOrder');
    minAx = min(cell2mat(Z'));
    maxAx = max(cell2mat(Z'));
    thetas = (0:N)*2*pi/N;
    for k = plotRange
        % tracker
        subplot(1,2,1);
        cla; grid on; hold on; axis equal;
        title(sprintf('time step %d', k));
        gated = tracker.gate(Z{k}, probbar(:, k), xbar(:, :, k), Pbar(:, :, :, k));
        minG = 1e20 * ones(2,1);
        maxG = zeros(2,1);
        [skupd, xkupd, Pkupd] = tracker.conditionalUpdate(Z{k}(:, gated), probbar(:, k), xbar(:, :, k), Pbar(:, :, :, k));
        beta = tracker.associationProbabilities(Z{k}(:, gated), probbar(:, k), xbar(:, :, k), Pbar(:, :, :, k));
        for s = 1:2
            plot(squeeze(xhat(1, s, 1:(k-1))), squeeze(xhat(2, s, 1:(k-1)))','Color', co(s,:));
            axis([minAx(1), maxAx(1), minAx(2), maxAx(2)])
            for j = 1:size(xkupd, 3)
                csj = mTL * co(s, :) + (1 - mTL) * (beta(j)*skupd(s, j)*co(s, :) + (1 - beta(j)*skupd(s, j)) * ones(1, 3)); % transparancy
                pl = plot([squeeze(xhat(1, s, k-1)); xkupd(1, s, j)], [squeeze(xhat(2, s, k-1)); xkupd(2, s, j)], '--', 'Color', csj);
                axis([minAx(1), maxAx(1), minAx(2), maxAx(2)])
                %alpha(pl, beta(j)*skupd(s, j));
                drawnow;
                covEllData = chol(Pkupd(1:2,1:2, s, j))' * [cos(thetas); sin(thetas)] + xkupd(1:2, s, j);
                pl = plot(covEllData(1, :), covEllData(2, :), '--', 'Color', csj);
            end

            [vk, Sk] = imm.modeFilters{s}.innovation([0,0], squeeze(xbar(:, s, k)), squeeze(Pbar(:, :, s, k)));
            gateData = chol(Sk)' * [cos(thetas); sin(thetas)] * sqrt(tracker.gateSize) + squeeze(xbar(1:2, s, k));
            plot(gateData(1, :),gateData(2, :), '.--', 'Color', co(s,:))
            scatter(Z{k}(1, :), Z{k}(2, :), 'rx')
            if a(k) > 0
                scatter(Z{k}(1, a(k)), Z{k}(2, a(k)), 'gx')
            end
            minGs = min(gateData, [], 2);
            minG = minGs .* (minGs < minG) + minG .* (minG < minGs);
            maxGs = max(gateData, [], 2);
            maxG = maxGs .* (maxGs > maxG) + maxG .* (maxG > maxGs);

        end
        scale = 1;
        minAx = minG - scale * (maxG - minG);
        maxAx = maxG + scale * (maxG - minG);
        axis([minAx(1), maxAx(1), minAx(2), maxAx(2)])
        %legend()

        % mode probabilities
        subplot(1,2,2)
        cla; grid on; hold on;
        for s = 1:2
            plot(probhat(s, 1:(k - 1)), 'Color', co(s, :));
            for j = 1:size(xkupd, 3)
                csj = mTL * co(s, :) + (1 - mTL) * (beta(j)*skupd(s, j)*co(s, :) + (1 - beta(j)*skupd(s, j)) * ones(1, 3)); % transparancy
                plot([k-1, k], [probhat(s, k-1), skupd(s, j)], '--', 'color', csj)
            end
        end
        axis([1, plotRange(end), 0, 1])
        drawnow;
        pause(plotpause)
    end
end

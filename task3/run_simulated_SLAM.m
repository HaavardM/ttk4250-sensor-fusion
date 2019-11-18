load simulatedSLAM;
K = numel(z);
%%
Q = ...
R = ...
doAsso = true;
JCBBalphas = [..., ...] % first is for joint compatibility, second is individual 
slam = EKFSLAM(Q, R, doAsso, JCBBalphas);

% allocate
xpred = cell(1, K);
Ppred = cell(1, K);
xhat = cell(1, K);
Phat = cell(1, K);
a = cell(1, K);

% init
xpred{1} = poseGT(:,1); % we start at the correct position for reference
Ppred{1} = zeros(3, 3); % we also say that we are 100% sure about that


figure(10); clf;
axAsso = gca;
N = K;
doAssoPlot = true; % set to true to se the associations that are done
for k = 1:N
    [xhat{k}, Phat{k}, NIS(k), a{k}] =  slam.update(xpred{k}, Ppred{k}, z{k});
    if k < K
        [xpred{k + 1}, Ppred{k + 1}] = slam.predict(xhat{k}, Phat{k}, odometry(:, k));
    end
    
    
    % checks
    if size(xhat{k},1) ~= size(Phat{k},1)
        error('dimensions of mean and covariance do not match')
    end
    
    if doAssoPlot && k > 1 %&& any(a{k} == 0) % uncoment last part to only see new creations
        cla(axAsso); hold on;grid  on;
        zpred = reshape(slam.h(xpred{k}), 2, []);
        scatter(axAsso, z{k}(1, :), z{k}(2, :));
        scatter(axAsso, zpred(1, :), zpred(2, :));
        plot(axAsso, [z{k}(1, a{k}>0); zpred(1, a{k}(a{k}>0))], [z{k}(2, a{k}>0); zpred(2, a{k}(a{k}>0))], 'r', 'linewidth', 2)
        
        legend(axAsso, 'z', 'zbar', 'a')
        title(axAsso, sprintf('k = %d: %s', k, sprintf('%d, ',a{k})));
        pause();
    end
end

% plotting
figure(3);
k = N;
clf;
%subplot(1,2,1);
hold on;

scatter(landmarks(1,:), landmarks(2,:), 'r^')
scatter(xhat{k}(4:2:end), xhat{k}(5:2:end), 'b.')

lh1 = plot(poseGT(1, 1:k), poseGT(2,1:k), 'r', 'DisplayName', 'gt');
lh2 = plot(cellfun(@(x) x(1), xhat), cellfun(@(x) x(2), xhat), 'b', 'DisplayName', 'est');

el = ellipse(xhat{k}(1:2),Phat{k}(1:2,1:2),5,200);
plot(el(1,:),el(2,:),'b');

for ii=1:((size(Phat{k}, 1)-3)/2)
   rI = squeeze(Phat{k}(3+[1,2]+(ii-1)*2,3+[1,2]+(ii-1)*2));
   el = ellipse(xhat{k}(3 + (1:2) + (ii-1)*2),rI,5,200);
   plot(el(1,:),el(2,:),'b');
end

axis equal;
title('results')
legend([lh1, lh2])
grid on;

% subplot(1,2,2);
% hold on;
% % funF = @(delta) sum([poseGT(1:2,k) - rotmat2d(delta(3)) * xhat(1:2,k) - delta(1:2), landmarks - rotmat2d(delta(3)) * reshape(xhat(4:end,k),2,[]) - delta(1:2)].^2 ,[1,2]);
% % funS = @(delta) sum([poseGT(1:2,1) - delta(1:2), landmarks - rotmat2d(delta(3)) * reshape(xhat(4:end,1),2,[]) - delta(1:2)].^2 ,[1,2]);
% funF = @(delta) sum([poseGT(1:2,k) - rotmat2d(delta(3)) * xhat{k}(1:2) - delta(1:2), landmarks(:, a) - rotmat2d(delta(3)) * reshape(xhat{k}(4:end),2,[]) - delta(1:2)].^2 ,[1,2]);
% funS = @(delta) sum([poseGT(1:2,1) - delta(1:2), landmarks - rotmat2d(delta(3)) * reshape(xhat{1}(4:end),2,[]) - delta(1:2)].^2 ,[1,2]);
% % deltaXFinal = fminunc(funF, zeros(3,1))
% % deltaXStart = fminunc(funS, zeros(3,1))
% 
% Rot = rotmat2d(deltaXFinal(3));
% %xcomp = [Rot, zeros(2,1); zeros(1, 2), 1]*xhat(1:3, :) + deltaXFinal;
% for k = 1:K
%     xcomp(:, k) = [Rot, zeros(2,1); zeros(1, 2), 1]*xhat{k}(1:3) + deltaXFinal;
% end
% %mcomp = Rot * reshape(xhat(4:end,k),2,[]) + deltaXFinal(1:2);
% mcomp = Rot * reshape(xhat{k}(4:end),2,[]) + deltaXFinal(1:2);
% 
% scatter(landmarks(1,:), landmarks(2,:), 'r^')
% scatter(mcomp(1,:), mcomp(2,:), 'b.')
% plot(poseGT(1, 1:k), poseGT(2,1:k), 'r');
% plot(xcomp(1,:), xcomp(2,:), 'b');
% 
% %el = ellipse(xcomp(1:2,k), Rot * squeeze(Phat(1:2,1:2,k)) * Rot',5,200);
% el = ellipse(xcomp(1:2,k), Rot * squeeze(Phat{k}(1:2,1:2)) * Rot',5,200);
% plot(el(1,:),el(2,:),'b');
% 
% for ii=1:m
%    %rI = squeeze(Rot * Phat(3+[1,2]+(ii-1)*2,3+[1,2]+(ii-1)*2,k) * Rot'); 
%    rI = squeeze(Rot * Phat{k}(3+[1,2]+(ii-1)*2,3+[1,2]+(ii-1)*2) * Rot'); 
%    el = ellipse(mcomp(:,ii),rI,5,200); 
%    plot(el(1,:),el(2,:),'b');
% end
% axis equal;
% grid on;
% title(sprintf('transformed: x = %0.2fm, y = %0.2fm , \\theta = %0.2fdeg',deltaXFinal(1:2),deltaXFinal(3)*180/pi))
%

%% consistency: what to do with variable measurement size..?
alpha = 0.05;
ANIS = mean(NIS)
ACI = chi2inv([alpha/2; 1 - alpha/2], 1)/N % NOT CORRECT NOW
CI = chi2inv([alpha/2; 1 - alpha/2], 1); % NOT CORRECT NOW
warning('These consistency intervals have wrong degrees of freedom')

figure(5); clf;
hold on;
plot(1:N, NIS(1:N));
insideCI = mean((CI(1) < NIS) .* (NIS <= CI(2)))*100;
plot([1, N], (CI*ones(1, 2))','r--');

title(sprintf('NIS over time, with %0.1f%% inside %0.1f%% CI', insideCI, (1-alpha)*100));
grid on;
ylabel('NIS');
xlabel('timestep');

%% run a movie
pauseTime = 0.05;
fig = figure(4);
ax = gca;
for k = 1:N
    cla(ax); hold on;
    scatter(ax, landmarks(1,:), landmarks(2,:), 'r^')
    scatter(ax, xhat{k}(4:2:end), xhat{k}(5:2:end), 'b*')
    plot(ax, poseGT(1, 1:k), poseGT(2,1:k), 'r-o','markerindices',10:10:k);
    plot(ax, cellfun(@(x) x(1), xhat(1:k)), cellfun(@(x) x(2), xhat(1:k)), 'b-o','markerindices',10:10:k);
    
    if k > 1 % singular cov at k = 1
        el = ellipse(xhat{k}(1:2),Phat{k}(1:2,1:2),5,200);
        plot(ax,el(1,:),el(2,:),'b');
    end
    
    for ii=1:((size(Phat{k}, 1)-3)/2)
       rI = squeeze(Phat{k}(3+[1,2]+(ii-1)*2,3+[1,2]+(ii-1)*2)); 
       el = ellipse(xhat{k}(3 + (1:2) + (ii-1)*2),rI,5,200);
       plot(ax, el(1,:),el(2,:),'b');
    end
    
    title(ax, sprintf('k = %d',k))
    grid(ax, 'on');
    pause(pauseTime);
end
        
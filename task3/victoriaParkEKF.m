%% Clear
clear all; close all

%% Meta stuff
nis_plot_mask = [true, true]; % [show plot, print plot]
results_plot_mask = [true, true];
results_bg_plot_mask = [true, true];
results_clustered_plot_mask = [true, true];
covar_plot_mask = [true, true];

show_and_print_nis_colored_track = false; % must be shown to be printed
show_lmk_count_update_time_plot = false; % do not automatically print this
show_lsr_movie = false;
doPlot = true; % show live slam movie

%% data loading
addpath('./victoria_park');
load('aa3_dr.mat');
load('aa3_lsr2.mat');
load('aa3_gpsx.mat');
timeOdo = time/1000; clear time;
timeLsr = double(TLsr)/1000; clear TLsr;
timeGps = timeGps/1000;
K = numel(timeOdo);
mK = numel(timeLsr);

%% View laser measurement movie
if show_lsr_movie
    viewLsr();
end

%% Parameters
% the car parameters
car.L = 2.83; % axel distance
car.H = 0.76; % center to wheel encoder
car.a = 0.95; % laser distance in front of first axel
car.b = 0.5; % laser distance to the left of center

doAsso = true;
doLAdd = true; % Add new landmarks

% the SLAM parameters
sigmas = [4e-2 , 4e-2, 2e-2];
CorrCoeff = [1, 0, 0; 0, 1, 0.9; 0, 0.9, 1];
Q = diag(sigmas) * CorrCoeff * diag(sigmas); % (a bit at least) emprically found, feel free to change

R = diag([1e-2, 1e-2]);

JCBBalphas = [1e-5, 1e-3]; % first is for joint compatibility, second is individual 
sensorOffset = [car.a + car.L; car.b];

% initialize TWEAK THESE TO BETTER BE ABLE TO COMPARE TO GPS
init_angle = 37;
init_offset_x = -1;
init_offset_y = 0;
P = zeros(3,3); % we say that we start knowing where we are in our own local coordinates

N = 15000;

%% Slam dunk
% allocate all the stuff we need
eta(1:3) = [Lo_m(1) + init_offset_x; La_m(2) + init_offset_y; init_angle * pi /180]; % set the start to be relatable to GPS. 
xupd = zeros(3, mK);
a = cell(1, mK);
slam = EKFSLAM(Q, R, doAsso, doLAdd, JCBBalphas, sensorOffset);
mk = 2; % first seems to be a bit off in timing
t = timeOdo(1);
updatecomptimes = []; % contains computational time of update steps
updatelandmarkcount = []; % contains number of landmarks after each update

if doPlot
    % Live slam plot
    figure(1); clf;  hold on; grid on; axis equal;
    ax = gca;
    % cheeper to update plot data than to create new plot objects
    scatter(Lo_m(timeGps < timeOdo(N)), La_m(timeGps < timeOdo(N)), '.')
    shLmk = scatter(ax, nan, nan, 'rx');
    shZ = scatter(ax, nan, nan, 'bo');
    lhPose = plot(ax, eta(1), eta(2), 'k', 'LineWidth', 2);
    th = title(ax, 'start');

    % Show covariance image
    figure, handlecovarfig = axes;
end
for k = 1:N
    if mk < mK && timeLsr(mk) <= timeOdo(k+1)
        dt = timeLsr(mk) - t;
        if  dt >= 0
            t = timeLsr(mk);
            odo = odometry(speed(k + 1), steering(k + 1), dt, car);
            [eta, P] = slam.predict(eta, P, odo);
        else
            error('negative time increment...')
        end
        z = detectTreesI16(LASER(mk,:));
        tic;
        [eta, P, NIS(k), a{k}] = slam.update(eta, P, z);
        updatecomptimes = [updatecomptimes; toc];
        updatelandmarkcount = [updatelandmarkcount; (size(eta, 1) - 3) / 2];
        xupd(:, mk) = eta(1:3); 
        NISmk(mk) = NIS(k);
        mk = mk + 1;
        if doPlot
            lhPose.XData = [lhPose.XData, eta(1)];
            lhPose.YData = [lhPose.YData, eta(2)];
            shLmk.XData = eta(4:2:end);
            shLmk.YData = eta(5:2:end);
            if ~isempty(z)
                zinmap = rotmat2d(eta(3)) * (z(1,:).*[cos(z(2,:)); sin(z(2,:))] + slam.sensOffset) + eta(1:2);
                shZ.XData = zinmap(1,:);
                shZ.YData = zinmap(2,:);
            end

            th.String = sprintf('step %d, laser scan %d, landmarks %d, measurements %d, num new = %d', k, mk, (size(eta,1) - 3)/2, size(z,2), nnz(a{k} == 0));
            image(P/trace(Q), 'Parent', handlecovarfig);
            drawnow;
            pause(0.01)
        end
    end
    
    if k < K
        dt = timeOdo(k+1) - t;
        t = timeOdo(k+1);
        odo = odometry(speed(k+1), steering(k+1), dt, car);
        [eta, P] = slam.predict(eta, P, odo);
    end
    if mod(k, 50) == 0
        fprintf(1, "Completed %d/%d timesteps\n", k, N);
    end
end

%% Plot results
if results_plot_mask(1)
    fig = figure(2);
elseif results_plot_mask(2)
    fig = figure("visible", "off");
end
if any(results_plot_mask)
    clf;  hold on; grid on; axis equal;
    plot(xupd(1, 1:(mk-1)), xupd(2, 1:(mk-1)))
    scatter(Lo_m(timeGps < timeOdo(N)), La_m(timeGps < timeOdo(N)), '.')
    scatter(eta(4:2:end), eta(5:2:end), 'rx');
    xlim([-150, 150]);
    ylim([-150, 150]);
    if results_plot_mask(2)
        printplot(fig, "a3-real-results.pdf");
    end
end
%% Plot results with background
if any(results_bg_plot_mask)
    if results_bg_plot_mask(1)
        fig = figure(19);
    elseif resutls_bg_plot_mask(2)
        fig = figure("visible", "off");
    end
    clf;  hold on; grid on; axis equal;
    scatter(eta(4:2:end), eta(5:2:end), 'y.');
    scatter(Lo_m(timeGps < timeOdo(N)), La_m(timeGps < timeOdo(N)), 'r.')
    plot(xupd(1, 1:(mk-1)), xupd(2, 1:(mk-1)), 'b')
    xlim([-150, 150]);
    ylim([-150, 150]);
    I = imread('aerial-wide.png');
    I = imrotate(I, -3, 'bilinear');
    impos_x = [-540 400];
    impos_y = [315 -167];
    implot = image(I, 'XData', impos_x, 'YData', impos_y);
    uistack(implot, 'bottom');
    if results_bg_plot_mask(2)
        printplot(fig, "a3-real-results-bg.pdf");
    end
end
%% Plot results with clustered landmarks
if any(results_clustered_plot_mask)
    if results_clustered_plot_mask(1)
        fig = figure(42);
    elseif results_clustered_plot_mask(2)
        fig = figure("visible", "off");
    end
    clf;  hold on; grid on; axis equal;
    plot(xupd(1, 1:(mk-1)), xupd(2, 1:(mk-1)))
    scatter(Lo_m(timeGps < timeOdo(N)), La_m(timeGps < timeOdo(N)), '.')
    c1_start = 170;
    c1_end = 306;
    c2_start = 570;
    c2_end = 630;
    scatter(eta(4:2:end), eta(5:2:end), 'k.');
    scatter(eta(c1_start:2:c1_end), eta(c1_start+1:2:c1_end+1), 'g*');
    scatter(eta(c2_start:2:c2_end), eta(c2_start+1:2:c2_end+1), 'r*');
    xlim([-150, 150]);
    ylim([-150, 150]);
    if results_clustered_plot_mask(2)
        printplot(fig, "a3-real-results-clustered.pdf");
    end
end
%% Plot NIS
if any(nis_plot_mask)
    if nis_plot_mask(1)
        fig = figure(3);
    elseif nis_plot_mask(2)
        fig = figure("visible", "off");
    end
    clf;
    plot(NISmk); hold on;
    line([1, numel(NISmk)],[1,1], 'Color', 'red', 'LineStyle', '--')
    line([1, numel(NISmk)],[0,0], 'Color', 'red', 'LineStyle', '--')
    insideCI = round(100*mean((0 <= NISmk).* (NISmk <= 1))); % NIS scaled to be between 0 and 1
    title(sprintf("NIS divided by chi2 upper bound (%d percent inside 95-CI)", insideCI));
    if nis_plot_mask(2)
        printplot(fig, "a3-real-nis.pdf");
    end
end

%% Plot NIS-colored track
if show_and_print_nis_colored_track
    plotcoloredtrack(xupd(1:2, 1:(mk-1)) + 2, NISmk, "NIS colored track", 3, 4);
    printplot(gcf, "a3-nis-colored-track.pdf");
end

%% Show covariance matrix
if any(covar_plot_mask)
    if covar_plot_mask(1)
        fig = figure(5);
    elseif covar_plot_mask(2)
        fig = figure("visible", "off");
    end
    image(P/trace(Q));
    if covar_plot_mask(2)
        printplot(fig, "a3-real-covar-matrix.pdf");
    end
end

%% Plot landmark count vs update computational time
if show_lmk_count_update_time_plot
    figure(6); clf; hold on;
    myfit = fit(updatelandmarkcount,updatecomptimes,'poly2','Robust','on');
    scatter(updatelandmarkcount, updatecomptimes, 'b.');
    legend("Update times");
    fitplt = plot(myfit, 'r--');
    fitplt.LineWidth = 2;
    xlabel("Number of landmarks");
    ylabel("Computation time of update [s]");
    ylim([0, 1]);
end
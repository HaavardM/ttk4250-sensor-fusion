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

%% Parameters
% the car parameters
car.L = 2.83; % axel distance
car.H = 0.76; % center to wheel encoder
car.a = 0.95; % laser distance in front of first axel
car.b = 0.5; % laser distance to the left of center

% the SLAM parameters
sigmas = ...
CorrCoeff = [1, 0, 0; 0, 1, 0.9; 0, 0.9, 1];
Q = diag(sigmas) * [1, 0, 0; 0, 1, 0.9; 0, 0.9, 1] * diag(sigmas); % (a bit at least) emprically found, feel free to change

R = ...

JCBBalphas = [..., ...]; % first is for joint compatibility, second is individual 
sensorOffset = [car.a + car.L; car.b];
slam = EKFSLAM(Q, R, true, JCBBalphas, sensorOffset);

% allocate
xupd = zeros(3, mK);
a = cell(1, mK);

% initialize TWEAK THESE TO BETTER BE ABLE TO COMPARE TO GPS
eta = [Lo_m(1); La_m(2); 30 * pi /180]; % set the start to be relatable to GPS. 
P = zeros(3,3); % we say that we start knowing where we are in our own local coordinates

mk = 2; % first seems to be a bit off in timing
t = timeOdo(1);
tic
N = 15000;


doPlot = true;
figure(1); clf;  hold on; grid on; axis equal;
ax = gca;
% cheeper to update plot data than to create new plot objects
lhPose = plot(ax, eta(1), eta(2), 'k');
shLmk = scatter(ax, nan, nan, 'rx');
shZ = scatter(ax, nan, nan, 'bo');
th = title(ax, 'start');

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
        [eta, P, NIS(k), a{k}] = slam.update(eta, P, z);
        xupd(:, mk) = eta(1:3); 
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
        toc
        disp(k)
        tic
    end
end

%% 
figure(2); clf;  hold on; grid on; axis equal;
plot(xupd(1, 1:(mk-1)), xupd(2, 1:(mk-1)))
scatter(Lo_m(timeGps < timeOdo(N)), La_m(timeGps < timeOdo(N)), '.')
scatter(eta(4:2:end), eta(5:2:end), 'rx');

% what can we do with consistency..? divide by the number of associated
% measurements? 
figure(3); clf;
plot(NIS);

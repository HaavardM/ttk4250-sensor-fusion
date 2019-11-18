function odo = odometry(ve, alpha, dt, car)
%% Odometry model for the victoria park dataset
% uses the  ackerman model for converting the center steering angle alpha
% and the left rear wheel encoder to speed at the center of the rear axel
%   $v_c = v_e/(1 - H\tan(\alpha)/L$
% where H is the distance from the center to the wheel, and L is the
% distance between the axels.
% 
% the dynamic model is written as 
%   $\dot{x} = v_c \cos(\phi)$
%   $\dot{y} = v_c \sin(\phi)$
%   $\dot{\phi} = \omega = v_c/r = v_c \tan(\alpha)/L$

% ve: the speed given in the data set
% alpha: the steering given in the dataset
vc = ve./(1 - car.H * tan(alpha)/car.L);
dp = dt * vc .* tan(alpha) /car.L;
dx = dt * vc .* sinc(dp/pi);
if abs(dp) < 0.001
    dy = dt * vc .* (dp/2 - dp.^3/24 + dp.^5/720);
else
    dy = dt * vc .* (cos(dp) - 1)./dp;
end

odo = [dx; dy; dp];
end

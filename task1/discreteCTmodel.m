function model = discreteCTmodel(q, r)
    % return a model structure for CT with position measurements that
    % implements the dynamic and measurement function along with their
    % Jacobian, as well as the process and measurement noise.
    % q;    q(1) is acceleration noise covariance.
    %       q(2) is turn rate noise covariance-
    % r;    positional measurement covariance
    model.f = @ f_m2_withT; %f_CT;
    model.F = @ Phi_m2_withT; %F_CT;
    model.Q = @ (x, Ts) [q(1) * [   Ts^3/3, 0,      Ts^2/2,	0;
                    0,      Ts^3/3, 0,    	Ts^2/2 ;
                    Ts^2/2, 0,      Ts,   	0,     ;
                    0,      Ts^2/2, 0,    	Ts,    ],   zeros(4,1);
                    zeros(1, 4),                        q(2) * Ts];
                
   model.h = @(x) x(1:2);
   model.H = @(x) [eye(2), zeros(2,3)];
   model.R = @(x) r*eye(2);
end

function x = f_CT(x,T)
    if(abs(x(5)) > 0.0001)
        x = [x(1) + sin(T * x(5)) * x(3) / x(5) - (1 - cos(T * x(5))) * x(4) / x(5);
                x(2) + (1 - cos(T * x(5))) * x(3) / x(5) + sin(T * x(5)) * x(4) / x(5);
                cos(T * x(5)) * x(3) - sin(T * x(5)) * x(4);
                sin(T * x(5)) * x(3) + cos(T * x(5)) * x(4);...
                x(5)];
    else
        x = [x(1) + T*x(3); x(2) + T*x(4); x(3); x(4); 0];
    end
end

function F = F_CT(x, T)
    if (abs(x(5)) > 0.0001)
        omega = x(5);
        sTo = sin(T * omega);
        cTo = cos(T * omega);
        F5 = [  (T * cTo * omega - sTo)/omega^2,     (-T * sTo * omega + 1 - cTo)/omega^2;
                (T * sTo * omega - 1 + cTo)/omega^2, (T * cTo * omega - sTo)/omega^2;
                -T * sTo,                            -T*cTo;
                T * cTo,                             -T*sTo;] * x(3:4);
             
        F = [[1, 0, sTo/omega, -(1-cTo)/omega;
              0, 1, (1-cTo)/omega, sTo/omega;
              0, 0, cTo, -sTo;
              0, 0, sTo, cTo], F5;
              zeros(1, 4), 1];
    else
        F = [eye(4) + diag(T*ones(2,1), 2), [-T^2 * x(4)/2; T^2 * x(3)/2; -T*x(4); T*x(3)];
            zeros(1,4), 1];
    end
end

function xout = f_m2_withT(x,T)

%Should now have been readjusted for parametrisation used in PDA

if(abs(x(5)) > 0.0001)
    xout = [x(1) + sin(T*x(5))*x(3)/x(5)-(1-cos(T*x(5)))*x(4)/x(5);...
        x(2) + (1-cos(T*x(5)))*x(3)/x(5)+sin(T*x(5))*x(4)/x(5);...
        cos(T*x(5))*x(3) - sin(T*x(5))*x(4);...
        sin(T*x(5))*x(3) + cos(T*x(5))*x(4);...
        x(5)];
else
    xout = [x(1) + T*x(3); x(2) + T*x(4); x(3); x(4); 0];
end
end

function Linmatrix = Phi_m2_withT(x,T)

if(abs(x(5)) > 0.0001)
        
    Jacobi_omega = [cos(T*x(5))*T*x(3)/x(5) - sin(T*x(5))*x(3)/x(5)^2 ...
        - sin(T*x(5))*T*x(4)/x(5) + (1-cos(T*x(5)))*x(4)/x(5)^2;...               %x    
        sin(T*x(5))*T*x(3)/x(5) - (1-cos(T*x(5)))*T*x(3)/x(5) ...
        + cos(T*x(5))*T*x(4)/x(5) - sin(T*x(5))*x(4)/x(5)^2;...                     %y         
        - sin(T*x(5))*T*x(3) - cos(T*x(5))*T*x(4);...                               %v_x
        cos(T*x(5))*T*x(3) - sin(T*x(5))*T*x(4);...                                 %v_y
        1];
    
    r= x(5);
    %u = x(3);
    %v = x(4);
    
    colX = [1,0,0,0,0]';
    colY = [0,1,0,0,0]';
    colU = [sin(r*T)/r,(1-cos(r*T))/r,cos(r*T),sin(r*T),0]';
    colV = [-(1-cos(r*T))/r,sin(r*T)/r,-sin(r*T),cos(r*T),0]';
    
    Linmatrix = [colX,colY,colU,colV, Jacobi_omega];
    
    
    
%     Linmatrix = cat(2,[1, sin(T*x(5))/x(5), 0, - (1-cos(T*x(5)))/x(5); ...
%         0, cos(T*x(5))*T, 0, -sin(T*x(5))*T; ...    
%         0, (1-cos(T*x(5)))*T/x(5), 1, sin(T*x(5))*T/x(5); ...
%         0, sin(T*x(5))*T, 0, cos(T*x(5))*T; ...
%         0, 0, 0, 0] , Jacobi_omega(x))
    
    %B = [1,0,0,0,0;0,0,1,0,0;0,1,0,0,0;0,0,0,1,0;0,0,0,0,1];
    %Linmatrix = B*Linmatrix*B';
else
    
    Linmatrix = [1, 0, T, 0, -T^2*x(4)/2; ...
        0, 1, 0, T, T^2*x(3)/2; ...
        0, 0, 1, 0, -T*x(4); ...
        0, 0, 0, 1, T*x(3); ...
        0, 0, 0, 0, 1];
    
    
%     
%     Linmatrix = [1, T, 0, 0, -T^2*x(4)/2; ...
%         0, 1, 0, 0, -T*x(4); ...
%         0, 0, 1, T, T^2*x(3)/2; ...
%         0,0,0,1,T*x(3); ...
%         0, 0, 0, 0, 1];
end
end
function [X, Z, a] = simulate_atc_track(Ts, K, q, r, init, PD, lambda, theCase)

if theCase
    rng('default');
    rng(8); % 8 good
end

% some limits
maxSpeed = 20;
maxTurn = pi/4;

models = cell(2,1);
models{1}.f = @(x, Ts) [x(1:4);0] + [Ts * x(3:4); zeros(3,1)];
models{1}.Q = @(x, Ts) [q(1) * [Ts^3/3, 0,      Ts^2/2,     0;
                        0,      Ts^3/3, 0,          Ts^2/2;
                        Ts^2/2,  0,     Ts,        0;
                        0,      Ts^2/2, 0,         Ts], zeros(4,1);
                        zeros(1, 4), 1e-20];
models{1}.sqrtQ = @(x, Ts) chol(models{1}.Q(x, Ts))';


models{2}.f = @(x, Ts) f_m2_withT(x, Ts);
models{2}.Q = @(x, Ts) [q(1) * [   Ts^3/3, 0,      Ts^2/2,	0;
                    0,      Ts^3/3, 0,    	Ts^2/2 ;
                    Ts^2/2, 0,      Ts,   	0,     ;
                    0,      Ts^2/2, 0,    	Ts,    ],   zeros(4,1);
                    zeros(1, 4),                        q(2) * Ts];
models{2}.sqrtQ = @(x, Ts) chol(models{2}.Q(x, Ts))';

 
h = @(x) x(1:2);
R = r * eye(2);
sqrtR = chol(R)';
n = 5;

xzero = init.x + chol(init.P)' * randn(n, 1);

X = zeros(size(xzero,1),K);
m = size(R,1);
z = zeros(m,floor(K));

X(:,1) = xzero;
z(:,1) = h(X(:,1)) + sqrtR*randn(m,1);

s = 1;

for k=2:K
    maxSpeedRatio = (X(3, k)^2 + X(4, k)^2)/(maxSpeed^2);
    if maxSpeedRatio > 1
        X(3:4, k) = X(3:4, k) /sqrt(maxSpeedRatio);
    end
    % limit turn rate
    if abs(X(5, k)) > maxTurn
        X(5, k) = maxTurn * sign(X(5));
    end
    
    if k == 25  %First turn
        X(5, k - 1) = (pi/2)/(18 * Ts);
        s = 2;
    elseif k == 43 %Straight movement
        X(5, k - 1) = 0;
        s = 1;
    elseif k == 68 % second turn
        X(5, k - 1) = (-2*pi/3)/(15*Ts);
        s = 2;
    elseif k == 74 %Straight movement
        X(5, k - 1) = 0;
        s = 1;
    end
    
    X(:,k) = models{s}.f(X(:, k-1), Ts) + models{s}.sqrtQ(X(:, k - 1), Ts) * randn(n,1);
    z(:, k) = h(X(:,k))+ sqrtR*randn(m,1);
end


minZ = min(z, [], 2);
maxZ = max(z, [], 2);
minZ = minZ - 0.05*(maxZ - minZ);
maxZ = maxZ + 0.05*(maxZ - minZ);
A = prod(maxZ - minZ);

Z = cell(K, 1);
a = zeros(K, 1);
for k = 1:K
    phi_k = poissrnd(lambda * A);
    Zfalse = (maxZ - minZ) .* rand(2, phi_k) + minZ;
    if rand(1) < PD
        a(k) = randi(phi_k + 1, 1);
        Z{k} = [Zfalse(:, 1:(a(k) - 1)), z(:, k), Zfalse(:, a(k):end)];
    else
        Z{k} = Zfalse;
    end
end

end

function xout = f_m2_withT(x,T)
    if(abs(x(5)) > 0.0001)
        xout = [x(1) + sin(T * x(5)) * x(3) / x(5) - (1 - cos(T * x(5))) * x(4) / x(5);
                x(2) + (1 - cos(T * x(5))) * x(3) / x(5) + sin(T * x(5)) * x(4) / x(5);
                cos(T * x(5)) * x(3) - sin(T * x(5)) * x(4);
                sin(T * x(5)) * x(3) + cos(T * x(5)) * x(4);...
                x(5)];
    else
        xout = [x(1) + T*x(3); x(2) + T*x(4); x(3); x(4); 0];
    end
end
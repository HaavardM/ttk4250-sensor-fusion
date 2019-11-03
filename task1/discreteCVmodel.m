function model = discreteCVmodel(q, r)
    % returns a structure that implements a discrete time CV model with
    % continuous time accelleration covariance q and positional
    % measurement with noise with covariance r, both in two dimensions.
    model.f = @(x, Ts) [x(1:4); 0] + [Ts * x(3:4); zeros(3,1)];
    model.F = @(x, Ts) [eye(4) + diag(Ts * ones(2,1), 2), zeros(4,1);
                        zeros(1, 5)];
    model.Q = @(x, Ts) [q * [Ts^3/3, 0,      Ts^2/2,     0;
                            0,      Ts^3/3, 0,          Ts^2/2;
                            Ts^2/2,  0,     Ts,        0;
                            0,      Ts^2/2, 0,         Ts], zeros(4,1);
                            zeros(1, 5)];
                    
    model.h = @(x) x(1:2);
    model.H = @(x) [eye(2), zeros(2,3)];
    model.R = @(x) r * eye(2);
end


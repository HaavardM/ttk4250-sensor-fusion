classdef EKFSLAM
    properties
        Q
        R
        doAsso
        alpha
        sensOffset
    end
    methods
        function obj = EKFSLAM(Q, R, doAsso, alphas, sensorOffset)
            obj.Q = Q;
            obj.R = R;
            
            if nargin < 3
                doAsso = false;
            end
            obj.doAsso = doAsso;
            
            if nargin < 4
                alphas = [0.001, 0.0001];
            end
            obj.alpha = alphas;
            
            if nargin < 5
                sensorOffset = zeros(2,1);
            end
            obj.sensOffset = sensorOffset(:);
        end
        
        function xpred = f(~, x, u)
            xpred = %
        end
        
        function Fx = Fx(~, x, u)
            Fx = %
            
            % check that jacobian is correct, remove for speed
            if norm(F - jacobianFD(@(X) obj.f(X, u), x, 1e-5), 'fro') > 1e-3
                error('some error in pred Jac')
            end
        end
        
        function Fu = Fu(~, x, u)
              Fu = %
              
              % check that the jacobian is correct, remove for speed
            if norm(F - jacobianFD(@(U) obj.f(x, U), x, 1e-5), 'fro') > 1e-3
                error('some error in pred Jac')
            end            
        end
        
        function [etapred, P] =  predict(obj, eta, P, zOdo)
            x = eta(1:3); % pose
            m = eta(4:end); % map
            
            xpred = %
            Fx = %
            Fu = %
            
            % in place for performance
            P(1:3, 1:3) = %
            P(1:3, 4:end) = % 
            P(4:end, 1:3) = % 
            
            % concatenate pose and landmarks again
            etapred = [xpred; m];
            
            % check that the covariance makes sense
            if any(eig(Ppred) <= 0) % costly, remove when tested
                warn('EKFpredict got cov not PSD')
            end
        end
        
        function zpred = h(obj, eta)
            x = eta(1:3); % pose 
            m = reshape(eta(4:end), 2 ,[]); % map (2 x m now)
            
            Rot = rotmat2d(-x(3)); % rot from world to body
            
            % cartesian measurement in world
            z_c = ... - Rot' * obj.sensOffset;
            
            % in body
            z_b = ...
            
            % polar
            zpred = ...
            
            % make column again
            zpred = zpred(:); 
        end
        
        function H = H(obj, eta)
            x = eta(1:3); % pose
            m = reshape(eta(4:end), 2 ,[]); % map
            
            numM = size(m, 2);  % number of landmarks
            
            Rot = rotmat2d(x(3));
            
            m_minus_rho = ...
            z_c = ... - Rot * obj.sensOffset;
            
            zpred = reshape(obj.h(eta), 2, []);
            zr = zpred(1,:);
            
            Rpihalf = [0, -1; 1, 0]; 
            
            % allocate
            Hx  = zeros(2 * numM, 3); % pose columns
            Hm = zeros(2 * numM, 2 * numM); % map columns (the rest)
            
            for i = 1:numM
                inds = 2*(i - 1) + [1; 2];
                
                jac_z_b = ... 
                
                Hx(inds(1), :) = ... jac z_r 
                Hx(inds(2), :) = ... jac z_phi
            
                Hm(inds, inds) = ... should be negative of the two first colums of Hx
            end
        
            % concatenate the H matrix
            H = [Hx, Hm];
            
            % check that it is done correctly, remove for speed
            if norm(H - jacobianFD(@(X) obj.h(X), eta, 1e-5), 'fro') > 1e-3
                error('some error in meas Jac')
            end
        end
        
        function [etaadded, Padded] = addLandmarks(obj, eta, P, z)
            n = size(P, 1);
            numLmk = numel(z)/2;
            
            % allocate
            lmnew = zeros(size(z));
            Gx = zeros(numLmk * 2, 3);
            Rall = zeros(numLmk * 2, numLmk * 2);

            for j = 1:numLmk
                % find indeces and the relevant measurement
                inds = 2 * (j - 1) + [1, 2];
                zj = z(inds);
                
                rot = rotmat2d(zj(2) + eta(3));

                lmnew(inds) = ... % mean
                Gx(inds, :) = ... % jac h^1 wrt. x
                Gz =  % jac h^-1 wrt. z
                Rall(inds, inds) = ... % the linearized measurement noise

            end
            
            % augment state
            etaadded = ...
            
            % add covariances
            Padded = blkdiag(P, ...);
            Padded((n+1):end, 1:n) = ...
            Padded(1:n, (n+1):end) = ...

            % sanity check, remove for speed
            if any(eig(Pupd) <= 0) % costly, remove when tested
                warning('EKFupdate got cov not PSD after adding a landmark');
            end
        end
        
        function [z, zpred, H, S, a] = associate(obj, z, zpred, H, S)
            if obj.doAsso
                % associate
                a = JCBB(z, zpred, S, obj.alpha(1), obj.alpha(2));

                % extract associated measurements
                zinds = false(size(z));
                zinds(1:2:end) = a > 0;
                zinds(2:2:end) = zinds(1:2:end);
                z = z(zinds);

                % extract and rearange predicted measurements and cov
                zbarinds = reshape([2*a(a>0) - 1, 2*a(a>0)]', [], 1);
                zpred = zpred(zbarinds);
                S = S(zbarinds, zbarinds);
                H = H(zbarinds, :);
            % else  % if no association is to be done, assume that all measurements are there and in order
            end
        end
        
        function [etaupd, Pupd, NIS, a] = update(obj, eta, P, z)    
            numLmk = (numel(eta) - 3)/2; % number of landmarks
            if numLmk > 0
                % prediction and innovation covariance
                zpred = ...
                H = ...
                S = ...
                z = z(:); % vectorize
                
                % perform data association if it is asked for
                [za, zpred, H, S, a] = obj.associate(z, zpred, H, S);

                % create the associated innovation
                v = za(:) - zpred;
                v(2:2:end) = wrapToPi(v(2:2:end)); % angles are in [-pi, pi]

                % Kalman update
                W = ...
                etaupd = ...
                NIS = ...
                Pupd = ...
                
                % sanity check, remove for speed
                if any(eig(Pupd) <= 0) % costly, remove when tested
                    warn('EKFupdate got cov not PSD');
                end
            else % all measurements are new landmarks
                a = zeros(size(z, 2), 1);
                z = z(:);
                NIS = 0;
                etaupd = eta;
                Pupd = P;
            end
            
            % create new landmarks if any is available
            if obj.doAsso
                isNewLmk = (a == 0);
                if any(isNewLmk)
                    % extract unassociated measurements
                    zNewInds = false(size(z));
                    zNewInds(1:2:end) = isNewLmk;
                    zNewInds(2:2:end) = isNewLmk;
                    znew = z(zNewInds);
                    
                    % create new landmarks
                    [etaupd, Pupd] = ...
                end
            end  
        end
    end
end

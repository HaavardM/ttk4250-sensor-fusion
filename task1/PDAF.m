classdef PDAF
    properties
        ekf
        
        clutterRate % expected number of false alarms
        PD % detection probability
        gateSize
    end
    methods
        function obj = PDAF(stateFilter, clutterRate, PD, gateSize)
            % init function
            obj = obj.setModel(stateFilter, clutterRate, PD, gateSize);
        end
        
        function obj = setModel(obj, ekf, clutterRate, PD, gateSize)
            % sets the parameters of the class object
            obj.ekf = ekf;
            obj.clutterRate = clutterRate;
            obj.PD = PD;
            obj.gateSize = gateSize;
        end
        
        function [xp, Pp] = predict(obj, x, P, Ts)
            % predict state distribution
            %
            % x (n x 1): mean to predict
            % P (n x n): covariance to predict
            % Ts (scalar): sampling time
            %
            % xp (n x 1): predicted mean
            % Pp (n x n): predicted covariance
            
            [xp, Pp] = obj.ekf.predict(x, P, Ts);
        end
        
        function gated = gate(obj, Z, x, P)
            % gates/validates measurements: (z-h(x))'S^(-1)(z-h(x)) <= g^2
            %
            % Z (dim(z) x m): measurements to gate
            % x (n x 1): state mean
            % P (n x n): state covariance
            %
            % gated (m x 1): gated(j) = true if measurement j is within gate
            
            m = size(Z, 2);
            gated = false(m, 1);
            gSquared = obj.gateSize;
            for j = 1:m
                [v, S] = obj.ekf.innovation(Z(:, j), x, P);
                g = (v' * (S \ v)) ;
                gated(j) = g < gSquared;
            end
        end
        
        function ll = loglikelihoodRatios(obj, Z, x, P)
            % Calculates the poseterior event loglikelihood ratios.
            %
            % Z (dim(z) x m): measurements to use in likelihoods
            % x (n x 1): state mean
            % P (n x n): state covariance
            % 
            % ll (m + 1 x 1): the posterior log likelihood ratios, ll(1)
            %                 corresponds to no detection
            
            % precalculate some parameters
            m = size(Z, 2);
            logPD = log(obj.PD);
            logPND = log(1 - obj.PD); % P_ND = 1 - P_D
            logClutter = log(obj.clutterRate);
            
            % allocate
            llCond = zeros(m, 1); % log(l^a),
            ll = zeros(m + 1, 1);
            
            
            % calculate log likelihood ratios
            ll(1) = logPND + logClutter; 
            for j = 1:m
                llCond(j) = obj.ekf.loglikelihood(Z(:, j), x, P);
                ll(j + 1) =  logPD + llCond(j);
            end
        end
        
        function beta = associationProbabilities(obj, Z, x, P)
            % calculates the poseterior event/association probabilities
            % 
            % Z (dim(z) x m): measurements ot use to get probabilities
            % x (n x 1): state mean
            % P (n x n): state covariance
            % 
            % beta (m + 1 x 1): the association probabilities (normalized
            %                   likelihood ratios)
           
           %log likelihoods
           lls = obj.loglikelihoodRatios(Z, x, P);
           
           % probabilities
           beta = exp(lls);
           beta = beta / sum(beta); % should sum to 1
        end
        
        function [xupd, Pupd] = conditionalUpdate(obj, Z, x, P)
            % updates the state with all possible measurement associations
            %
            % Z (dim(z) x m): measurements to use for update
            % x (n x 1): state mean to update
            % P (n x n): state covariance to update
            %
            % xupd (n x m + 1): the updated states for all association
            %                   events. xupd(:, 1) corresponds to no detection
            % Pupd (n x n x m + 1): the updated covariances for all association
            %                   events. Pupd(:, :, 1) corresponds to no detection
           
            m = size(Z, 2);
            
            % allocate
            xupd = zeros([size(x, 1), m + 1]);
            Pupd = zeros([size(P), m + 1]);
            
            % undetected
            xupd(:, 1) = x;
            Pupd(:, :, 1) = P;
            
            % detected
            for j = 1:m 
               [xupd(:, j + 1), Pupd(:, :, j + 1)] = obj.ekf.update(Z(:, j), x, P);
            end
        end
        
        function [xred, Pred] = reduceMixture(obj, beta, x, P)
            % reduces a Gaussian mixture to a single Gauss
            % 
            % beta (m + 1 x 1): the mixture weights 
            % x (n x m + 1): the means to reduce
            % P (n x n x m + 1): the covariances to reduce
            %
            % xred (n x 1): the mean of the mixture
            % Pred (n x n): the covariance of the mixture
            
            [xred, Pred] = reduceGaussMix(beta, x, P); %... % Hint: reduceGaussMix from assignment 3
        end
        
        function [xupd, Pupd] = update(obj, Z, x, P)
            % The whole PDAF update sycle.
            %
            % Z (dim(z) x m): measurements to use for update
            % x (n x 1): state mean to update
            % P (n x n): state covariance to update
            %
            % xupd (n x 1): the mean of the PDAF update
            % Pupd (n x n): the covariance of the PDAF update
            
            % remove the not gated measurements from consideration
            gated = obj.gate(Z, x, P); % ... 
            Zg = Z(:, gated);
            
            % find association probabilities
            beta = obj.associationProbabilities(Zg, x, P); % ...
            
            % find the mixture components pdfs
            [xcu, Pcu] = obj.conditionalUpdate(Zg, x, P); %...
            
            % reduce mixture
            [xupd, Pupd] = obj.reduceMixture(beta, xcu, Pcu); %...
        end
    end
end
    
function lse = logSumExp(a)
    % more numerically stable way(less chance of underflow and overflow)
    % to calculate logsumexp of the list a.
    % 
    % uses the fact
    % log(sum(exp(a))) = log(sum(exp(b)exp(a - b))
    % = log(exp(b)sum(exp(a - b))) = b + log(sum(exp(a - b)))
    % where we let b = max(a), 
    amax = max(a(:));
    lse = amax + log(sum(exp(a - amax)));
end
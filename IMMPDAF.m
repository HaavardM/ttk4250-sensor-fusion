classdef IMMPDAF
    properties
        imm
        
        clutterRate % expected number of false alarms
        PD % detection probability
        gateSize
    end
    methods
        function obj = IMMPDAF(imm, clutterRate, PD, gateSize)
            % init function
            obj = obj.setModel(imm, clutterRate, PD, gateSize);
        end
        
        function obj = setModel(obj, imm, clutterRate, PD, gateSize)
            % sets the parameters of the class object
            obj.imm = imm;
            obj.clutterRate = clutterRate;
            obj.PD = PD;
            obj.gateSize = gateSize;
        end
        
        function [sprobsp, xp, Pp] = predict(obj, sprobs, x, P, Ts)
            % predict state distribution
            %
            % sprobs (M x 1): mode probabilities to predict
            % x (n x M): means to predict
            % P (n x n x M): covariances to predict
            % Ts (scalar): sampling time
            %
            % sprobsp (M x 1): predicted mode probabilities
            % xp (n x M): predicted means
            % Pp (n x n x M): predicted covariances
            
            [sprobsp, xp, Pp] = %... 
        end
        
        function gated = gate(obj, Z, sprobs, x, P)
            % gates/validates measurements per mode: (z-h(x))'S^(-1)(z-h(x)) <= g^2
            %
            % Z (dim(z) x m): measurements to gate
            % sprobs (M x 1): mode probabilities
            % x (n x M): state mean
            % P (n x n x M): state covariances
            %
            % gated (m x 1): gated(j) = true if measurement j is within
            %       gate of at least one mode
            
            m = size(Z, 2);
            gated = false(m, 1);
            gSquared = obj.gateSize;
            for j = 1:m
                [NIS, NISes] = obj.imm.NIS(Z(:,j), sprobs, x, P); 
                gated(j) = %...
            end
        end
        
        function ll = loglikelihoodRatios(obj, Z, sprobs, x, P)
            % Calculates the posterior event loglikelihood ratios.
            %
            % Z (dim(z) x m): measurements to use in likelihoods
            % sprobs (M x 1): mode probabilities
            % x (n x M): state mean
            % P (n x n x M): state covariance
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
            ll(1) = % association loglikelihood ratio for no detection
            for j = 1:m
                llCond(j) = %... calculate imm loglikelihood
                ll(j + 1) = %... association loglikelihood ratio for detection j
            end
        end
        
        function beta = associationProbabilities(obj, Z, sprobs, x, P)
            % calculates the poseterior event/association probabilities
            % 
            % Z (dim(z) x m): measurements ot use to get probabilities
            % sprobs (M x 1): mode probabilities
            % x (n x M): state means
            % P (n x n x M): state covariances
            % 
            % beta (m + 1 x 1): the association probabilities (normalized
            %                   likelihood ratios)
           
           %log likelihoods
           lls = obj.loglikelihoodRatios(Z, sprobs, x, P);
           
           % probabilities
           beta = %... 
        end
        
        function [sprobsupd, xupd, Pupd] = conditionalUpdate(obj, Z, sprobs, x, P)
            % updates the state with all possible measurement associations
            %
            % Z (dim(z) x m): measurements to use for update
            % sprobs (M x 1): mode probabilities to update
            % x (n x M): state means to update
            % P (n x n x M): state covariances to update
            %
            % sprobsupd (M x m + 1): the updated mode probabilities for all
            %                   association events. (:, 1) correspongs to no detection
            % xupd (n x M x m + 1): the updated states for all association
            %                   events. xupd(:, :, 1) corresponds to no detection
            % Pupd (n x n x M x m + 1): the updated covariances for all association
            %                   events. Pupd(:, :, :, 1) corresponds to no detection
           
            m = size(Z, 2);
            
            % allocate
            sprobsupd = zeros([size(sprobs, 1), m + 1]);
            xupd = zeros([size(x), m + 1]);
            Pupd = zeros([size(P), m + 1]);
            
            % undetected
            sprobsupd(:, 1) = %... 
            xupd(:, :, 1) = %...
            Pupd(:, :, :, 1) = %...
            
            % detected
            for j = 1:m 
               [sprobsupd(:, j + 1), xupd(:, :, j + 1), Pupd(:, :, :, j + 1)] = %... update conditioned on measurement j
            end
        end
        
        function [sprobsred, xred, Pred] = reduceMixture(obj, beta, sprobs, x, P)
            % reduces a Gaussian mixture to a single Gauss
            % 
            % beta (m + 1 x 1): the association probabilities (mixture weights)
            % sprobs (M x m + 1): the mixture conditional mode probabilities
            % x (n x M x m + 1): the means to reduce
            % P (n x n x M x m + 1): the covariances to reduce
            %
            % sprobsred (M x 1): the marginal mode probabilities 
            % xred (n x M): the mean of the mode mixtures
            % Pred (n x n x M): the covariance of the mode mixtures
            
            M = size(sprobs, 1);
            
            joint = %.. Joint probability for mode and association (M x m + 1)
            sprobsred = %... marginal mode probabilities (M x 1)
            betaCondS = %... association probabilites conditionend on the mode probabilites (M x m + 1)
            
            xSize = size(x);
            PSize = size(P);
            xred = zeros(xSize(1:2));
            Pred = zeros(PSize(1:3));
            for s = 1:M
                [xred(:, s), Pred(:, : ,s)] = %... mean and variance per mode
            end
        end
        
        function [sprobsupd, xupd, Pupd] = update(obj, Z, sprobs, x, P)
            % The whole PDAF update sycle.
            %
            % Z (dim(z) x m): measurements to use for update
            % sprobs (M x 1): mode probabilities to update
            % x (n x 1): state mean to update
            % P (n x n): state covariance to update
            %
            % sprobsupd (M x 1): the updated mode probabilities
            % xupd (n x 1): the mean of the PDAF update
            % Pupd (n x n): the covariance of the PDAF update
            
            % remove the not gated measurements from consideration
            gated = %...
            Zg = Z(:, gated);
            
            % find association probabilities
            beta = %...
            
            % find the mixture components (conditional update)
            [sprobscu, xcu, Pcu] = %...
            
            % reduce mixture
            [sprobsupd, xupd, Pupd] = %...
        end
    end
end
    
function lse = logSumExp(a)
    % more numerically stable way(less chance of underflow and overflow)
    % to calculate logsumexp of a list, a.
    % 
    % uses the fact
    % log(sum(exp(a))) = log(sum(exp(b)exp(a - b))
    % = log(exp(b)sum(exp(a - b))) = b + log(sum(exp(a - b)))
    % where we let b = max(a), 
    amax = max(a(:));
    lse = amax + log(sum(exp(a - amax)));
end
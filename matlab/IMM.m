classdef IMM
   properties
       modeFilters  % cell of EKFs
       PI           % markov transition matrix.
       M            % number of modes
   end
   methods
       function obj = IMM(modelcellarr, PI)
           % modelcell (M x 1 cell): cell array of EKFs 
           % PI (M x M): Markov transition matrix
            obj = obj.setModel(modelcellarr, PI);
       end
       
       function obj = setModel(obj, modelcellarr, PI)
               % sets the internal functions and paramters
               %
               % modelcell (M x 1 cell): cell array of EKFs 
               % PI (M x M): Markov transition matrix
               obj.modeFilters = modelcellarr;
               obj.PI = PI;
               obj.M = size(PI, 1);
       end
       
       
       function [spredprobs, smixprobs] = mixProbabilities(obj, sprobs)
           % IMM: step 1
           %
           % probs (M x 1): mode probabilities
           %
           % spredprobs (M x 1): predicted mode probabilities
           % smixprobs (M x M): mixing probabilities

           % Joint probability for this model and next
           spsjointprobs = obj.PI .* (ones(obj.M, 1) * sprobs(:)');

           % marginal probability for next model
           spredprobs = sum(spsjointprobs, 2);

           % conditionional probability for model at this time step on the next.
           smixprobs = spsjointprobs ./ (spredprobs * ones(1, obj.M));
       end

       function [xmix, Pmix] = mixStates(obj, smixprobs, x, P)
           % IMM: step 2
           % smixprob (M x M): mixing probabilities
           % x (dim(state) x M): means to mix
           % P (dim(state) x dim(state) x M): covariances to mix
           %
           % xmix (dim(state) x M): mixed means
           % Pmix (dim(state) x dim(state) x M): mixed covariances
           
           % allocate
           xmix = zeros(size(x));
           Pmix = zeros(size(P));
           
           % mix for each mode, 
           for i = 1:obj.M
               [xmix(:, i), Pmix(:, :, i)] = reduceGaussMix(smixprobs(i, :), x, P);
           end
       end
       
       function [xpred, Ppred] = modeMatchedPrediction(obj, x, P, Ts)
           % IMM: prediction part of step 3 
           % x (dim(state) x M matrix): mean to predict
           % P (dim(state) x dim(state) x M): covariance to predict
           % Ts: sampling time for prediction.
           %
           % xpred (dim(state) x M): predicted means
           % Ppred (dim(state) x dim(state) x M): predicted covariances
           
           % allocate
           xpred = zeros(size(x));
           Ppred = zeros(size(P));
           
           % mode matched prediction
           for i = 1:obj.M
               [xpred(:, i), Ppred(:, :, i)] = obj.modeFilters{i}.predict(x(:, i), P(:, :, i), Ts);
           end
       end
       
       function [sprobspred, xpred, Ppred] = predict(obj, sprobs, x, P, Ts)
           % IMM: step 1, 2 and prediction part of 3
           % sprobs (M x 1): mode probabilities
           % x (dim(state) x M): means to predict
           % P (dim(state) x dim(state) x M): covariances to predict
           % Ts: sampling time
           % 
           % sprobspred (M x 1): predicted mode probabilities
           % xpred (dim(state) x M): predicted means
           % Ppred (dim(state) x dim(state) x M): predicted covariances
           
           % step 1
           [sprobspred, smixprobs] = obj.mixProbabilities(sprobs);
           
           % step 2
           [xmix, Pmix] = obj.mixStates(smixprobs, x, P); % ...
           
           % prediction part of step 3
           [xpred, Ppred] = obj.modeMatchedPrediction(xmix, Pmix, Ts);
       end
       
       function [xupd, Pupd, logLambdas] = modeMatchedUpdate(obj, z, x, P)
           % IMM: update part of step 3
           % z (dim(measurement) x 1): measurement
           % x (dim(state) x M): the means to update
           % P (dim(state) x dim(state) x M): covariances to update
           %
           % xupd (dim(state) x M): updated means
           % Pupd (dim(state) x dim(state) x M): updated covariances
           % logLambdas (M x 1): measurement loglikelihood for given modes
           
           % allocate
           xupd = zeros(size(x));
           Pupd = zeros(size(P));
           logLambdas = zeros(obj.M, 1);
           
           % mode matched update and likelihood
           for i = 1:obj.M
               [xupd(:, i), Pupd(:, :, i)] = obj.modeFilters{i}.update(z, x(:, i), P(:, :, i));
               logLambdas(i) = obj.modeFilters{i}.loglikelihood(z, x(:, i), P(:, :, i));
           end
           
       end
       
       function [supdprobs, loglikelihood] = updateProbabilities(obj, logLambdas, sprobs)
           % IMM: step 4
           %
           % logLambdas (M x 1): measurement loglikelihood for given modes
           % sprobs (M x 1): mode probabilities
           %
           % supdprobs (M x 1): updated mode probabilities
           % loglikelihood: measurement log likelilhood (total, ie. p(z_k | z_(1:k-1)))
           
           ... % you might want to do some precalculations here.
           
           loglikelihood = logSumExp(logLambdas + log(sprobs)); % you might want to use the logSumExp function at the bottom of this file
           supdprobs = exp(logLambdas + log(sprobs) - loglikelihood);
       end
       
       function [supdprobs, xupd, Pupd, loglikelihood] = update(obj, z, sprobs, x, P)
           % IMM: combining update part of step 3 and step 4
           %
           % z (dim(measurement) x 1): measurement
           % sprobs (M x 1): mode probabilities
           % x (dim(state) x M): the means to update
           % P (dim(state) x dim(state) x M): covariances to update
           %
           % supdprobs (M x 1): updated mode probabilities
           % xupd (dim(state) x M): updated means
           % Pupd (dim(state) x dim(state) x M): updated covariances
           % loglikelihood: measurement log likelilhood (total, ie. p(z_k | z_(1:k-1)))
           
           % update part of step 3
           [xupd, Pupd, logLambdas] = obj.modeMatchedUpdate(z, x, P);
           
           % step 4
           [supdprobs, loglikelihood] = obj.updateProbabilities(logLambdas, sprobs);
       end
       
       function [xest, Pest] = estimate(obj, sprobs, x, P)
           % IMM: step 5. A single mean and covariance as estimate. Reuse
           % of reduceGaussMix should simplify things.
           %
           % sprobs (M x 1): mode probabilities
           % x (dim(state) x M): means per mode 
           % P (dim(state) x dim(state) x M): covariances per mode
           %
           % xest (dim(state) x M): MMSE/mean estimate
           % Pest (dim(state) x dim(state) x M): covariance of the estimation error.
           
           [xest, Pest] = reduceGaussMix(sprobs, x, P);
       end
       
       function [NIS, NISes] = NIS(obj, z, sprobs, x, P)
           % calculate the NIS for each mode, and one for the averaged
           % innvoations.
           %
           % sprobs (M x 1): mode probabilities
           % x (dim(state) x M): means per mode 
           % P (dim(state) x dim(state) x M): covariances per mode
           %
           % NIS (scalar): NIS calculated based on the estimation mean and covariance
           % NISes (M x 1): NIS for each mode
           
           m = size(z,1);
           NISes = zeros(obj.M, 1);
           innovs = zeros(m, obj.M);
           Ss = zeros(m, m, obj.M);
           for s = 1:obj.M
               [innovs(:, s), Ss(:, :, s)] = obj.modeFilters{s}.innovation(z, x(:, s), P(:, :, s));
               NISes(s) = obj.modeFilters{s}.NIS(z, x(:, s), P(:,:,s));
           end
           [totInnov, totS] = reduceGaussMix(sprobs, innovs, Ss);
           NIS = totInnov' * (totS \ totInnov);
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
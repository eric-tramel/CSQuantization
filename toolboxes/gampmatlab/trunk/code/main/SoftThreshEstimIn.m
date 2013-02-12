classdef SoftThreshEstimIn < EstimIn
    % SoftThreshEstimIn:  Inputs a soft thresholding scalar input function.
    % Allows GAMP to be used to solve "min_x 1/2/var*norm(y-A*x,2)^2 + lambda*norm(x,1)". 
    
    properties
        %lambda = the gain on the ell1 term in the MAP cost. 
	%The soft threshold is set according to the expression thresh = lambda * mur;
        lambda;
        maxSumVal = true;   % Max-sum GAMP (true) or sum-product GAMP (false)?
    end
    
    methods
        % Constructor
        function obj = SoftThreshEstimIn(lambda, maxSumVal)
            obj = obj@EstimIn;
            obj.lambda = lambda;
            if nargin >= 2 && ~isempty(maxSumVal) && isscalar(maxSumVal)
                obj.maxSumVal = logical(maxSumVal);
            end
        end
        
        % Prior mean and variance
        function [mean0, var0, valInit] = estimInit(~)
	    %Set these to reasonable constants
            mean0 = 0; 
            var0  = 5e-4;
            valInit = -inf;
        end
        
        % Carry out soft thresholding
        function [xhat,xvar,val] = estim(obj,rhat,rvar)
            if ~obj.maxSumVal
                % Compute sum-product GAMP updates
                
                % Begin by computing various constants on which the
                % posterior mean and variance depend
                sig = sqrt(rvar);                           % Gaussian std dev
                muL = rhat + obj.lambda.*rvar;          	% Lower integral mean
                muU = rhat - obj.lambda.*rvar;          	% Upper integral mean
                muL_over_sig = muL ./ sig;
                muU_over_sig = muU ./ sig;
                cdfL = normcdf(-muL_over_sig);              % Lower cdf
                cdfU = normcdf(muU_over_sig);               % Upper cdf
                NormConL = obj.lambda/2 .* ...             % Mass of lower integral
                    exp( (muL.^2 - rhat.^2) ./ (2*rvar) ) .* cdfL;
                NormConU = obj.lambda/2 .* ...             % Mass of upper integral
                    exp( (muU.^2 - rhat.^2) ./ (2*rvar) ) .* cdfU;
                NormCon = NormConL + NormConU;      % Posterior normaliz. constant recip.
                
                % Compute the ratio normpdf(a)/normcdf(a) for
                % appropriate upper- and lower-integral constants, a
                RatioL = 2/sqrt(2*pi) ./ erfcx(muL_over_sig / sqrt(2));
                RatioU = 2/sqrt(2*pi) ./ erfcx(-muU_over_sig / sqrt(2));
                
                % Now compute the first posterior moment...
                xhat = (1./NormCon) .* (NormConL .* ...
                    (muL - sig.*RatioL) + NormConU .* ...
                    (muU + sig.*RatioU));
                
                % ...and second central posterior moment
                varL = rvar .* (1 - RatioL.*(RatioL - muL./sig));
                varU = rvar .* (1 - RatioU.*(RatioU + muU./sig));
                meanL = muL - sig.*RatioL;
                meanU = muU + sig.*RatioU;
                SecondMoment = (1./NormCon) .* ...
                    ( NormConL .* (varL + meanL.^2) + ...
                    NormConU .* (varU + meanU.^2) );
                xvar = SecondMoment - xhat.^2;
                
                % Lastly, compute negative KL divergence:
                % \int_x p(x|r) log(p(x)/p(x|r)) dx
                val = (1/2)*log(2*pi*rvar) + ...
                    log( NormConL + NormConU ) + ...
                    (1./(2*rvar)).*(xvar + xhat.^2) - ...
                    rhat.*xhat./rvar + (rhat.^2)./(2*rvar);
            else
                % Compute max-sum GAMP updates
                
                %Compute the thresh
                thresh = obj.lambda .* rvar;
                
                %Estimate the signal
                xhat = max(0,abs(rhat)-thresh) .* sign(rhat);
                
                %Estimate the variance
                %xvar = rvar .* (mean(double(abs(xhat) > 0))*ones(size(xhat)));
                xvar = rvar .* (abs(xhat) > 0);
                
                %Output negative cost
                %val = -1*obj.lambda*abs(rhat);
                val = -1*obj.lambda.*abs(xhat);	% seems to work better
            end
        end
        
        % Computes p(y) for y = x + v, with x ~ p(x), v ~ N(0,yvar)
        function py = plikey(obj, y, yvar)
            mu = y;
            sig2 = yvar;
            sig = sqrt(sig2);                           % Gaussian prod std dev
            muL = mu + obj.lambda.*sig2;                % Lower integral mean
            muU = mu - obj.lambda.*sig2;                % Upper integral mean
            muL_over_sig = muL ./ sig;
            muU_over_sig = muU ./ sig;
            cdfL = normcdf(-muL_over_sig);              % Lower cdf
            cdfU = normcdf(muU_over_sig);               % Upper cdf
            NormConL = obj.lambda/2 .* ...              % Mass of lower integral
                exp( (muL.^2 - mu.^2) ./ (2*sig2) ) .* cdfL;
            NormConU = obj.lambda/2 .* ...              % Mass of upper integral
                exp( (muU.^2 - mu.^2) ./ (2*sig2) ) .* cdfU;
            py = NormConL + NormConU;       % Posterior normaliz. constant recip. 
        end
        
    end
    
end


classdef GaussMixEstimOut < EstimOut
    % GaussMixEstimOut:  Gaussian Mixture output estimation
    %
    % Corresponds to an output channel of the form
    %   y = z + q
    %
    % Where q has the density
    % p(q) = (1-lambda) Nor(0,nu0) + lambda Nor(0,nu1)
    
    properties
        
        % Measured data
        Y;
        
        %Variances
        nu0;
        nu1;
        
        %on probability
        lambda;
        

        
        
    end
    
    methods
        % Constructor
        function obj = GaussMixEstimOut(Y,nu0,nu1,lambda)
            obj = obj@EstimOut;
            obj.Y = Y;
            obj.nu0 = nu0;
            obj.nu1 = nu1;
            obj.lambda = lambda;

            
            %Warn user that logLike has not been implemented yet
            warning(['The logLike method is not yet implemented for this' ...
                ' class. It should not be used with adaptive step size'...
                ' until this problem is rectified']) %#ok<WNTAG>
        end
        
        
        
        % Estimation function
        % Computes the posterior mean and variance given the estimates phat
        % and pvar
        function [zhat, zvar] = estim(obj, phat, pvar)
            
            
            %Compute the intrinsic LLR
            int_llr = log(obj.lambda./(1 - obj.lambda));
            
            %Compute the extrinsic LLR
            ext_llr = (1/2)*log( (obj.nu0 + pvar) ./ (obj.nu1 + pvar)) ...
                + (obj.Y - phat).^2 .*...
                ( 1 ./ sqrt(obj.nu0 + pvar) - 1./ sqrt(obj.nu1 + pvar));
            
            %Limit ext_llr
            ext_llr = min(10,max(-10,ext_llr));
            
            %Now compute p1 and p0
            p1 = 1 ./ (1 + exp(-int_llr - ext_llr) + eps);
           
            
            %Compute p0
            p0 = 1 - p1;
            
            %We can now obtain the mean
            zeta0 = (obj.Y .* pvar + phat.*obj.nu0) ./ ...
                (obj.nu0 + pvar);
            
            zeta1 = (obj.Y .* pvar + phat.*obj.nu1) ./ ...
                (obj.nu1 + pvar);
            
            zhat = p0.*zeta0 + p1.*zeta1;
            
            %Compute variance
            zvar = p0.*(obj.nu0.*pvar ./ (obj.nu0 + pvar) + zeta0.^2)...
                + p1.*(obj.nu1.*pvar ./ (obj.nu1 + pvar) + zeta1.^2)...
                - zhat.^2;
            
            
            %Protect form zvar getting too small
            rval = 0.99;
            zvar = max(min(zvar,rval*pvar),-rval*pvar);
            
            
            
            
        end
        
        
        
        % Compute log likelihood
        %   E( log p_{Y|Z}(y|z) )
        function ll = logLike(~,zhat,~)
            
            %Not implemented yet
            ll = zeros(size(zhat));
            
            
        end
        
        function S = numColumns(obj)
            %Return number of columns of Y
            S = size(obj.Y,2);
        end
    end
    
end


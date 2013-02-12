classdef CAwgnEstimOut < EstimOut
    % CAwgnEstimOut:  CAWGN scalar output estimation function
    %
    % Corresponds to an output channel of the form
    %   y = scale*z + CN(0, wvar)
    
    properties
        % Prior mean and variance
        y;      % Measured output
        wvar;   % Variance
        scale = 1;  % scale factor
        
        % True indicates to compute output for max-sum
        maxSumVal = false;
    end
    
    methods
        % Constructor
        function obj = CAwgnEstimOut(y, wvar, maxSumVal, scale)
            obj = obj@EstimOut;
            obj.y = y;
            obj.wvar = wvar;
            if (nargin >= 3)
                if (~isempty(maxSumVal))
                    obj.maxSumVal = maxSumVal;
                end
            end
            if (nargin >= 4)
                obj.scale = scale;
            end
            
            %Warn user if they set zero variance
            if any(obj.wvar == 0)
                warning(['Tiny non-zero variances will be used for'...
                    ' computing log likelihoods. May cause problems'...
                    ' with adaptive step size if used.']) %#ok<*WNTAG>
            end
        end
        
        % AWGN estimation function
        % Provides the posterior mean and variance of variable z
        % from an observation y = scale*z + w, z = CN(zmean0,zvar0), w = CN(0,wvar)
        function [zmean, zvar] = estim(obj, zmean0, zvar0)
            
            % Compute posterior mean and variance
            s = obj.scale;
            gain = conj(s)*zvar0./((abs(s)^2)*zvar0 + obj.wvar);
            zmean = gain.*(obj.y-s*zmean0) + zmean0;
            zvar = obj.wvar.*zvar0./((abs(s)^2)*zvar0 + obj.wvar);
            
        end
        
        % Compute log likelihood
        % For max-sum GAMP, compute
        %   E( log p_{Y|Z}(y|z) ) with z = CN(zhat, zvar)
        % For sum-product compute
        %   log p_{Y|Z}(y|z) @ z = zhat
        function ll = logLike(obj,zhat,zvar)
            
            % Ensure variance is small positive number
            wvar1 = max(eps, obj.wvar);
            
            % Get scale
            s = obj.scale;            

            % Compute log-likelihood
            if ~(obj.maxSumVal)
                predErr = (abs(obj.y-s*zhat).^2 + (abs(s)^2)*zvar)./wvar1;
            else
                predErr = (abs(obj.y-s*zhat).^2)./wvar1;
            end
            ll = -(predErr); %return the values without summing
        end
        
        function S = numColumns(obj)
            %Return number of columns of Y
            S = size(obj.y,2);
        end
    end
    
end


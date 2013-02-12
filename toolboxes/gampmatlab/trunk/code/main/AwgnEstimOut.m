classdef AwgnEstimOut < EstimOut
    % AwgnEstimOut:  AWGN scalar output estimation function
    %
    % Corresponds to an output channel of the form
    %   y = scale*z + N(0, wvar)
    
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
        function obj = AwgnEstimOut(y, wvar, maxSumVal, scale)
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
		if (scale <= 0),
                  error('Fourth argument of AwgnEstimOut must be positive');
		end;
            end
            
            %Warn user about inputs
            if any(~isreal(obj.y)),
              error('First argument of AwgnEstimOut must be real-valued');
            end;
            if any((obj.wvar<0))||any(~isreal(obj.wvar)),
              error('Second argument of AwgnEstimOut must be non-negative');
            end;
            if any(obj.wvar==0)
                warning(['Tiny non-zero variances will be used for'...
                    ' computing log likelihoods. May cause problems'...
                    ' with adaptive step size if used.']) %#ok<*WNTAG>
            end
        end
        
        % AWGN estimation function
        % Provides the posterior mean and variance of variable z
        % from an observation y = scale*z + w, z = N(zmean0,zvar0), w = N(0,wvar)
        function [zmean, zvar] = estim(obj, zmean0, zvar0)
            
            % Compute posterior mean and variance
            s = obj.scale;
            gain = conj(s)*zvar0./((abs(s)^2)*zvar0 + obj.wvar);
            zmean = gain.*(obj.y-s*zmean0) + zmean0;
            zvar = obj.wvar.*zvar0./((abs(s)^2)*zvar0 + obj.wvar);
            
        end
        
        % Compute log likelihood
        % For max-sum GAMP, compute
        %   E( log p_{Y|Z}(y|z) ) with z = N(zhat, zvar)
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
            ll = -0.5*(predErr); %return the values without summing
        end
        
        function S = numColumns(obj)
            %Return number of columns of Y
            S = size(obj.y,2);
        end
    end
    
end


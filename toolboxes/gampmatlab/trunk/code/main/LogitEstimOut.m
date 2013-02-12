% CLASS: LogitEstimOut
% 
% HIERARCHY (Enumeration of the various super- and subclasses)
%   Superclasses: EstimOut
%   Subclasses: N/A
% 
% TYPE (Abstract or Concrete)
%   Concrete
%
% DESCRIPTION (High-level overview of the class)
%   The LogitEstimOut class defines a scalar observation channel, p(y|z),
%   that constitutes a logit binary classification model, i.e., y is an
%   element of the set {0,1}, z is a real number, and
%           Pr{y=0|z} = 1 - Pr{y=1|z} = 1 / (1+exp(scale*z)),
%   where scale is a scalar parameter that determines the relative shape of
%   the sigmoidal curve.
%
%   LogitEstimOut can operate in either of two different modes
%   corresponding to two different flavors of GAMP: sum-product mode (for
%   MMSE estimation) or max-sum mode (for MAP estimation).  In both modes,
%   approximations are required for the purpose of tractability.  In
%   sum-product mode, an intractable integral is replaced by a numerical
%   integration.  The number of discrete points included in the integration
%   is determined by the parameter Npts, while the range of discrete points
%   is the interval (-Wmax, Wmax), where Wmax is a real-valued parameter.
%   In max-sum mode, the maximization of a concave function is approximated
%   by a maximization of a second-order Taylor approximation of that
%   function.  An advantage of the max-sum mode is that it is quicker than
%   the alternative sum-product version.
%
% PROPERTIES (State variables)
%   y           An M-by-1 array of binary ({0,1}) class labels for the
%               training data, where M is the number of training data
%               points
%   scale       Scalar controlling the shape of the sigmoidal function.
%               Higher values produce sharper sigmoids [Default: 1]
%   npts        Number of points used in discrete integration, for
%               sum-product GAMP [Default: 100]
%   wmax        Range of points included in discrete integration, [-wmax,
%               wmax], for sum-product GAMP [Default: 4]
%   maxSumVal   Perform MMSE estimation (false) or MAP estimation (true)?
%               [Default: false]
%
% METHODS (Subroutines/functions)
%   LogitEstimOut(y)
%       - Default constructor.  Assigns remaining properties to default
%         values
%   LogitEstimOut(y, scale)
%       - Optional constructor.  Sets both y and scale.
%   LogitEstimOut(y, scale, npts)
%       - Optional constructor.  Sets y, scale, and npts.
%   LogitEstimOut(y, scale, npts, wmax)
%       - Optional constructor.  Sets y, scale, npts, and wmax.
%   LogitEstimOut(y, scale, npts, wmax, maxSumVal)
%       - Optional constructor.  Sets y, scale, npts, wmax, and maxSumVal.
%   estim(obj, zhat, zvar)
%       - Provides the posterior mean and variance of a variable z when
%         p(y|z) is the logit model and maxSumVal = false (see 
%         DESCRIPTION), and p(z) = Normal(zhat,zvar).  When maxSumVal =
%         true, estim returns MAP estimates of each element of z, as well
%         as the second derivative of log p(y|z).
%

%
% Last change: 08/28/12
% Change summary: 
%       - Created (v0.1) (10/20/11; SR)
%       - Added numColumns method (08/16/12; JAZ)
%       - Added support for max-sum GAMP (v0.2) (08/28/12; JAZ)
% Version 0.2
%

classdef LogitEstimOut < EstimOut
    
    properties        
        y;                  % Vector of labels: 0 or 1
        scale = 1;          % scale factor
        npts = 100;         % number of points used in discrete integration
        wmax = 4;           % maximum value of integration
        maxSumVal = false; 	% Sum-product (false) or max-sum (true) GAMP
    end
    
    methods
        % *****************************************************************
        %                      CONSTRUCTOR METHOD
        % *****************************************************************
        function obj = LogitEstimOut(y, scale, npts, wmax, maxsumval)
            obj = obj@EstimOut;
            obj.y = y;
            if (nargin >= 2) && ~isempty(scale)
                obj.scale = scale;
            end
            if (nargin >= 3) && ~isempty(npts)
                obj.npts = npts;
            end
            if (nargin >= 4) && ~isempty(wmax)
                obj.wmax = wmax;
            end
            if (nargin >= 5) && ~isempty(maxsumval)
                if isscalar(maxsumval)
                    obj.maxSumVal = logical(maxsumval);
                else
                    error('maxSumVal must be a logical scalar')
                end
            end
        end
        
        
        % *****************************************************************
        %                           SET METHODS
        % *****************************************************************
        function obj = set.y(obj, y)
            if ~all((y(:) == 0) | (y(:) == 1))
                error('Elements of y must be binary (0,1)')
            else
                obj.y = y;
            end
        end
        
        
        % *****************************************************************
        %                          ESTIM METHOD
        % *****************************************************************
        % This function will compute the posterior mean and variance of a
        % random vector z with prior distribution N(zhat0, zvar0), given
        % observations y obtained through the separable channel model:
        % Pr{y(m)=0|z(m)} = 1 - Pr{y(m)=1|z(m)} = 1 / (1+exp(scale*z(m))),
        % if obj.maxSumVal = false, otherwise it will return zhat = argmax
        % log p(y|z) - 1/2/zvar0 (z - zhat0)^2 and the second derivative of
        % log p(y|z) evaluated at zhat in zvar, if obj.maxSumVal = true
        function [zhat, zvar] = estim(obj, zhat0, zvar0)
            
            % Check if zhat0 and zvar0 are only scalars (can occur during
            % first method call by gampEst) and resize
            if numel(zhat0) == 1
                zhat0 = zhat0*ones(size(obj.y));
                zvar0 = zvar0*ones(size(obj.y));
            end
            
            switch obj.maxSumVal
                case false
                    % Compute the sum-product GAMP updates
                    
                    % Gaussian pdf
                    w = linspace(-obj.wmax,obj.wmax,obj.npts)';
                    pw = exp(-0.5*w.^2);
                    
                    % Loop over points
                    ny = length(obj.y);
                    zhat = zeros(ny,1);
                    zvar = zeros(ny,1);
                    for iy = 1:ny
                        z = zhat0(iy) + sqrt(zvar0(iy))*w;
                        if (obj.y(iy))
                            pyz = 1./(1+exp(-obj.scale*z));
                        else
                            pyz = 1./(1+exp(obj.scale*z));
                        end
                        pzy = pyz.*pw;
                        pzy = pzy / sum(pzy);
                        zhat(iy) = z'*pzy;
                        zvar(iy) = ((z-zhat(iy)).^2)'*pzy;
                    end
                    
                case true
                    % Compute the max-sum GAMP updates
                    
                    PMonesY = 2*(obj.y - 1/2);  % Convert {0,1} to {-1,1}
                    a = obj.scale;              % Shorthand for scale
%                     NegPM = -PMonesY;           % Negated {-1,1}
                    
                    % Specify the expansion point for the Taylor series
                    % approximation
                    EP = (sign(PMonesY) == sign(zhat0)) .* zhat0;
                    
%                     % Find an approximation to 
%                     % zhat = argmax log p(y|z) - 1/2/zvar0 (z - zhat0)^2
%                     % by way of a second-order Taylor series approximation
% %                     numer = a * PMonesY .* zvar0 .* ...
% %                         (1 + exp(a * NegPM .* zhat0) );
% %                     denom = (a^2)*zvar0 + 2 + exp(-a*zhat0) + exp(a*zhat0);
%                     C = (1 + exp(a*PMonesY.*EP));
%                     numer = PMonesY .* ((EP - zhat0).*C - a*zvar0);
%                     denom = C + (a^2*zvar0.*exp(a*PMonesY.*EP) ./ C);
%                     zhat = EP - (numer ./ denom);
                    
                    % fsolve approach - Iteratively locate the value of z,
                    % zhat, that sets the derivative equal to zero
                    opts = optimset('Jacobian', 'on', 'MaxIter', 20, ...
                        'Display', 'off');
                    F = @(z) zero_deriv(obj, z, zhat0, zvar0);
                    zhat = fsolve(F, EP, opts);
                    
                    % Also compute the 2nd derivative (wrt z) of log p(y|z)
                    % evaluated at the zhat that was just computed
                    deriv = -a^2 ./ (2 + exp(a*zhat) + exp(-a*zhat));
                    
                    % Output in zvar a function of the 2nd derivative that,
                    % once manipulated by gampEst, yields the desired
                    % max-sum expression for -g'_{out}
                    zvar = zvar0 ./ (1 - zvar0.*deriv);
            end
       
        end
                
        
        % *****************************************************************
        %                         LOGLIKE METHOD
        % *****************************************************************
        % This function will compute *an approximation* to the expected
        % log-likelihood, E_z[log p(y|z)] when performing sum-product GAMP
        % (obj.maxSumVal = false).  The approximation is based on Jensen's 
        % inequality, i.e., computing log E_z[p(y|z)] instead.  If
        % performing max-sum GAMP (obj.maxSumVal = true), logLike returns
        % log p(y|z) evaluated at z = zhat
        function ll = logLike(obj,zhat,zvar)
            switch obj.maxSumVal
                case false
                    % Numerically integrate to estimate E_z[log p(y|z)]
                    
                    % Gaussian pdf
                    w = linspace(-obj.wmax,obj.wmax,obj.npts)';
                    pw = exp(-0.5*w.^2);
                    pw = pw / sum(pw);
                    
                    % Loop over points
                    ny = length(obj.y);
                    logpy = zeros(ny,1);
                    for iy = 1:ny
                        z = zhat(iy) + sqrt(zvar(iy))*w;
                        z = min(max(z,-10),10);
                        if (obj.y(iy))
                            logpy(iy) = -pw'*log(1+exp(-obj.scale*z));
                        else
                            logpy(iy) = -pw'*log(1+exp(obj.scale*z));
                        end
                    end
                    ll = sum(logpy);
                    
                case true
                    % Evaluate log p(y|zhat)
                    PMonesY = 2*(obj.y - 1/2);  % Convert {0,1} to {-1,1}
                    a = obj.scale;              % Shorthand for scale
                    NegPM = -PMonesY;           % Negated {-1,1}
                    ll = -log(1 + exp(a*NegPM.*zhat));
            end
        end
        
        
        % *****************************************************************
        %                       NUMCOLUMNS METHOD
        % *****************************************************************
        function S = numColumns(obj)
            %Return number of columns of y
            S = size(obj.y, 2);
        end
        
    end
    
    methods (Access = private)
        % *****************************************************************
        %                         FSOLVE METHOD
        % *****************************************************************
        % This method is used by MATLAB's fsolve function in the max-sum
        % case to locate the value of z that sets the prox-operator
        % derivative to zero
        function [Fval, Jacobian] = zero_deriv(obj, z, phat, pvar)
            % Compute value of derivative at z
            ExpZ = min(realmax, exp(-obj.scale*z));
            OhYeah = (1./pvar).*(z - phat);
            Fval = obj.scale*(obj.y - 1) - OhYeah + (obj.scale*obj.y - ...
                OhYeah).*ExpZ;
            
            % Optionally compute Jacobian of F at z
            if nargout >= 2
                Jvec = (obj.scale*OhYeah - (obj.scale)^2*obj.y - 1./pvar).*...
                    ExpZ - 1./pvar;
                Jacobian = diag(Jvec);      % Okay for small problems
            end
        end
    end
end


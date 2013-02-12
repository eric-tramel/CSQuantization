classdef EllpEstimIn < EstimIn
    % EllpEstimIn:  Inputs a soft thresholding scalar input function.
    % Allows GAMP to be used to solve "min_x 1/2/var*norm(y-A*x,2)^2 + lambda*norm(x,p)^p"
    % for 0<p<=1. 
    
    properties
        %lambda = the gain on the ellp term in the MAP cost. 
	%The soft threshold is set according to the expression thresh = lambda * mur;
        lambda;

	%Value of p defining the ell-p norm (0<p<=1, where p=1 for soft thresholding)
	p;
    end
    
    methods
        % Constructor
        function obj = EllpEstimIn(lambda,p)
            obj = obj@EstimIn;
            obj.lambda = lambda;
            obj.p = p;
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
            
            %Compute the thresh
            thresh = obj.lambda .* rvar;
            
            %Estimate the signal
            xhat = max(0,abs(rhat)-thresh.*abs(rhat).^(obj.p-1)) .* sign(rhat);
            
            %Estimate the variance
            xvar = rvar .* (1 - thresh.*(obj.p-1).*abs(rhat).^(obj.p-2)) .* (xhat~=0);
            
            %Output negative cost 
            %val = -1*obj.lambda*norm(rhat,obj.p)^obj.p;
%             val = -1*obj.lambda*norm(xhat,obj.p)^obj.p;		% seems to work better?
            val = -1*obj.lambda*abs(xhat).^obj.p;		% seems to work better?
            
        end
        
    end
    
end


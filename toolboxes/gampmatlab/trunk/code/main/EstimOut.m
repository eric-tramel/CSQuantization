classdef EstimOut < handle
% EstimOut:  Base class for output function in the G-AMP algorithm
%
% The output estimation function is defined based on the output channel
% probability distribution p_{Y|Z}(y|z).  The value y is not passed
% to the function, but typically stored as a member of the function.
    
    methods (Abstract)
        
        % Main estimation method:  For the sum-product algorithm,
        % the method should return:
        %
        %   zhat = E( Z | Y)
        %   zvar = var( Z | Y )
        %
        % where Z = N(phat, pvar)
        [zhat,zvar] = estim(obj,phat,pvar)
        
        % Log-likelihood:  The method should return 
        %   E( log p_{Y|Z}(y|Z) )  with Z = N(zhat,zvar)
        ll = logLike(obj,zhat,zvar) 
       
        
    end
    
    %Placeholder methods to be overwritten if desired
    methods
        
        %Return number of columns
        function S = numColumns(obj) %#ok<MANU>
            
            %By default, return 1
            S = 1;
        end
        
    end
    
end
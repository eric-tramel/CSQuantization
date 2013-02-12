classdef EstimInConcat < EstimIn
    % EstimInConcat:  Concatenation of input estimators
    %
    % The estimator creates an input estimator for a vector x = [x1 ... xn]
    % with estimators estimArray{i} 
    properties
        estimArray;  % Cell array of input estimators
        
        % Index array.  The components for xi are in x(ind(i):ind(i+1)-1)
        ind;         
    end
    
    methods
        
        % Constructor
        % estimArray is a cell array of estimators
        % nx(i) is the number of elements for the estimator estimArray{i}.
        function obj = EstimInConcat(estimArray, nx)    
            obj = obj@EstimIn;
            obj.estimArray = estimArray;
            
            % Constructor 
            nelem = length(estimArray);
            obj.ind = zeros(nelem+1,1);
            obj.ind(1) = 1;
            for ielem = 1:nelem
                obj.ind(ielem+1) = obj.ind(ielem)+nx(ielem);
            end
            
        end
                
        % Initial estimate 
        function [xhat, xvar, valInit] = estimInit(obj)
            
            % Initial estimates
            nelem = length(obj.estimArray);
            xhat = zeros(obj.ind(nelem+1)-1,1);
            xvar = zeros(obj.ind(nelem+1)-1,1);
            valInit = 0;
                        
            for i = 1:nelem
                
                % Get estimates
                I = (obj.ind(i):obj.ind(i+1)-1)';
                [xhati, xvari, valIniti] = obj.estimArray{i}.estimInit();
                
                % Store results
                xhat(I) = xhati;
                xvar(I) = xvari;
                valInit = valInit + sum(valIniti(:));                
            end
            
        end

        % Initial estimate based on the mean
        function [xhat, xvar, val] = estim(obj, rhat, rvar)
            
            % Initial estimates
            nelem = length(obj.estimArray);
            xhat = zeros(obj.ind(nelem+1)-1,1);
            xvar = zeros(obj.ind(nelem+1)-1,1);
            val = zeros(obj.ind(nelem+1)-1,1);            
                        
            for i = 1:nelem
                
                % Get estimates
                I = (obj.ind(i):obj.ind(i+1)-1)';
                [xhati, xvari, vali] = ...
                    obj.estimArray{i}.estim(rhat(I), rvar(I));
                
                % Store results
                xhat(I) = xhati;
                xvar(I) = xvari;
                val(I) = vali;                
            end
        end
        
    end
end


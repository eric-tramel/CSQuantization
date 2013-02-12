classdef EstimOutConcat < EstimOut
    % EstimOutConcat:  Concatenation of output estimators
    %
    % The estimator creates an output estimator for a vector z = [z1 ... zn]
    % with estimators estimArray{i} 
    properties
        estimArray;  % Cell array of input estimators
        
        % Index array.  The components for zi are in z(ind(i):ind(i+1)-1)
        ind;         
    end
    
    methods
        
        % Constructor
        % estimArray is a cell array of estimators
        % nz(i) is the number of elements for the estimator estimArray{i}.
        function obj = EstimOutConcat(estimArray, nz)    
            obj = obj@EstimOut;
            obj.estimArray = estimArray;
            
            % Constructor 
            nelem = length(estimArray);
            obj.ind = zeros(nelem+1,1);
            obj.ind(1) = 1;
            for ielem = 1:nelem
                obj.ind(ielem+1) = obj.ind(ielem)+nz(ielem);
            end
            
        end
                
        % Estimation function (see description in EstimOut.estim() )
        % Applied for each array element
        function [zhat,zvar] = estim(obj,phat,pvar)
            
            % Initial estimates
            nelem = length(obj.estimArray);
            zhat = zeros(obj.ind(nelem+1)-1,1);
            zvar = zeros(obj.ind(nelem+1)-1,1);
                        
            for i = 1:nelem
                
                % Get estimates
                I = (obj.ind(i):obj.ind(i+1)-1)';
                [zhati, zvari] = obj.estimArray{i}.estim(phat(I),pvar(I));
                
                % Store results
                zhat(I) = zhati;
                zvar(I) = zvari;
            end
            
        end

        % Log-likelihood applied to each component
        function ll = logLike(obj,zhat,zvar) 
            
            % Initial estimates
            nelem = length(obj.estimArray);
            ll = 0.0;
                        
            for i = 1:nelem
                
                % Get estimates
                I = (obj.ind(i):obj.ind(i+1)-1)';
                lli =  obj.estimArray{i}.logLike(zhat(I), zvar(I));
                
                % Accumulate log likelihood
                ll = ll + lli;
            end
        end
        
    end
end


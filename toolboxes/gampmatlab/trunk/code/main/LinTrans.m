classdef LinTrans < handle
    % LinTrans:  Abstract base class for specifying a linear transform in
    % the GAMP algorithm.  Linear transforms are specified by deriving from
    % this base class and implementing the methods listed below.  These
    % methods include operations to multiply by the matrix, its transpose
    % and the componentwise square of the matrix.
    %
    % The linear transform for the 
    methods (Abstract)
        
        % Size
        [m,n] = size(obj)
        
        % Matrix multiply:  z = A*x
        y = mult(obj,x)

        % Matrix multiply transpose:  x = A'*z
        x = multTr(obj,z)

        % Matrix multiply with square:  z = (A.^2)*x
        z = multSq(obj,x)

        % Matrix multiply with componentwise square transpose:  
        % x = (A.^2)'*z
        x = multSqTr(obj,z)
        
        %Optional method, not included in required interface. This method
        %returns an updated version of pvar based on matrix uncertainty in
        %the operator. gampEst checks using ismethod to see if this is
        %defined. Was not added to this template specifically to avoid
        %breaking existing code.
        %Include matrix uncertainty
        %function pvar = includeMatrixUncertainty(obj,pvarBar,xhat,xvar)
    end
end
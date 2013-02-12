% CLASS: GAMPState
% 
% HIERARCHY (Enumeration of the various super- and subclasses)
%   Superclasses: N/A
%   Subclasses: N/A
% 
% TYPE (Abstract or Concrete)
%   Concrete
%
% DESCRIPTION (High-level overview of the class)
%   This class is used to store information about the "state" of GAMP
%   after execution, which is required by certain functions of the TurboOpt
%   class in order to perform turbo inference.
%
% PROPERTIES (State variables)
%   * All properties of this class are hidden, and should not require user
%   access or manipulation *
%
% METHODS (Subroutines/functions)
%   * All methods of this class are hidden *
% 

%
% Coded by: Justin Ziniel, The Ohio State Univ.
% E-mail: zinielj@ece.osu.edu
% Last change: 12/09/11
% Change summary: 
%       - Created (12/08/11; JAZ)
% Version 1.0
%

classdef GAMPState

    properties (GetAccess = private, SetAccess = immutable, Hidden)
        xhat;
        xvar;
        rhat;
        rvar;
        shat;
        svar;
        zhat;
        zvar;
    end
    
    methods (Hidden)
        % Constructor
        function obj = GAMPState(xhat, xvar, rhat, rvar, shat, svar, ...
                zhat, zvar)
            if nargin < 8
                error('Insufficient number of input arguments')
            end
            obj.xhat = xhat;
            obj.xvar = xvar;
            obj.rhat = rhat;
            obj.rvar = rvar;
            obj.shat = shat;
            obj.svar = svar;
            obj.zhat = zhat;
            obj.zvar = zvar;
        end
        
        % Accessor
        function [xhat, xvar, rhat, rvar, shat, svar, zhat, zvar] = ...
                getState(obj)
            xhat = obj.xhat;
            xvar = obj.xvar;
            rhat = obj.rhat;
            rvar = obj.rvar;
            shat = obj.shat;
            svar = obj.svar;
            zhat = obj.zhat;
            zvar = obj.zvar;
        end
        
    end % methods
    
end % classdef
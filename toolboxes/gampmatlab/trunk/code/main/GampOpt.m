classdef GampOpt
    % Options for the GAMP optimizer.
    
    properties
        nit = 200;          % number of iterations
        
        % Remove mean from A by creating a new matrix with one additional
        % row and column.
        removeMean = false;
        
        % Relative minimum variance on output and input.
        % This prevents poor conditioning. The lower limit is relative to
        % the initial variance
        pvarMin = 1e-10;
        xvarMin = 1e-10;
        
        % Enable adaptive step size
        adaptStep = true;
        
        %Create a window for the adaptive step size test. Setting this to
        %zero causes it to have no effect. For postive integer values,
        %creates a moving window of this length when checking the step size
        %acceptance criteria. The new value is only required to be better
        %than the worst in this window, i.e. the step size criteria is not
        %required to monotonicaly increase. As with other algorithms, this
        %modification tends to improve convergence speed
        stepWindow = 20;
        
        %Set to true to use a Barzilai Borwein type rule to choose a new
        %step size after each succesful step. Otherwise, the step is set to
        %the previous succesful step size
        bbStep = false;
        
        %Option to use identical variances across coefficients. Allows 
	%alg to run about twice as fast, but may impact performance 
        uniformVariance = false;
        
        % Print progress
        verbose = false;
        
        %Specify a convergence tolerance. Iterations are terminated when
        %the norm of the differnece in two iterates divided by the norm of
        %the current iterate falls below this threshold. Set to -1 to
        %disable
        tol = 1e-4;
        
	% Initial step size
        step = 1;           

        % Minimum step size.  If step is initialized below this value,
	% then the iteration will be automatically considered successful.
        stepMin = 0;
        
        % Maximum allowed step size
        stepMax = 1;
        
        % Multiplicative step size increase, when successful
        stepIncr = 1;
        
        % Multiplicative step size decrease, when unsuccessful
        stepDecr = 0.5;
        
        % Iterations are termined when the step size becomes smaller
        % than this value. Set to -1 to disable
        stepTol = 1e-10;
        
        %Robust GAMP options
        %Logical flag to include a step size in the pvar/zvar calculation.
        %This momentum term often improves numerical performance. On by
        %defualt.
        pvarStep = true;
        
        %Option to "normalize" variances for computation. May improve
        %robustness in some cases, but fundamentally changes the quantities
	%that are slowed down when step<1
        varNorm = true;
        
        %Provide initial guesses for xhat0,xvar0,shat0. If these are set to
        %empty, then the appropriate method of the input estimator is
        %called to initialize these values. This functionality is useful
        %for warm starting the algorithm.
        xhat0 = [];
        xvar0 = [];
        shat0 = [];
        
	%The following are also useful for warm-starting when the
	%stepsize is not equal to one or varNorm=true.
	svar0 = [];
	xhatPrev0 = [];
	scaleFac = 1;
    end
    
    methods
        
        % Constructor with default options
        function opt = GampOpt(varargin)
            if nargin == 0
                % No custom parameters values, thus create default object
                return
            elseif mod(nargin, 2) == 0
                % User is providing property/value pairs
                names = fieldnames(opt);    % Get names of class properties
                
                % Iterate through every property/value pair, assigning 
                % user-specified values.  Note that the matching is NOT
                % case-sensitive
                for i = 1 : 2 : nargin - 1
                    if any(strcmpi(varargin{i}, names))
                        % User has specified a valid property
                        propName = names(strcmpi(varargin{i}, names));
                        opt.(propName{1}) = varargin{i+1};
                    else
                        % Unrecognized property name
                        error('GampOpt: %s is an unrecognized option', ...
                            num2str(varargin{i}));
                    end
                end
                return
            else
                error(['The GampOpt constructor requires arguments ' ...
                    'be provided in pairs, e.g., GampOpt(''adaptStep'',' ...
                    ' false, ''nit'', 50)'])
            end
        end
    end
    
end

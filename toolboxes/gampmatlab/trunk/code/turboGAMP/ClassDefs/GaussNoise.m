% CLASS: GaussNoise
% 
% HIERARCHY (Enumeration of the various super- and subclasses)
%   Superclasses: Observation
%   Subclasses: N/A
% 
% TYPE (Abstract or Concrete)
%   Concrete
%
% DESCRIPTION (High-level overview of the class)
%   This class is used to define the statistics of the additive white
%   Gaussian noise (AWGN) matrix, N, that corrupts the linear measurements
%                y(t) = A(t)*x(t) + w(t),   t = 1, ..., T
%   where y(t) is a real- or complex-valued length-M vector of 
%   measurements, and each element of w(t) is assumed to follow a real- or 
%   circularly-symmetric-complex-valued Gaussian distribution with zero 
%   mean, i.e.,
%               pdf(w(m,t)) = Normal(w(m,t); 0, PSI(m,t)),
%   for m = 1, ..., M, t = 1, ..., T.  Note that this class is also 
%   appropriate in the single measurement vector (SMV) case, (i.e., T = 1).
%
%   The distribution of W is specified by providing the noise variance, 
%   PSI.  Acceptable inputs include a scalar, an M-by-1 row vector, a 
%   1-by-T column vector, or an M-by-T matrix.  In the first three cases, 
%   the input is replicated to yield an M-by-T matrix of variances.  By
%   default, it is assumed that the noise is real-valued.  If the noise is
%   complex-valued, then set data = 'complex'.
%
%   To learn the noise variance from the data using an expectation-
%   maximization (EM) learning algorithm, the property learn_prior_var can
%   be set to 'scalar', to learn a single common noise variance for all
%   elements of W, 'row', to learn a unique noise variance for each row of
%   W, or 'column', to learn a unique noise variance for each column of W.
%   To prevent the EM algorithm from learning the noise variance, set
%   learn_prior_var to 'false'
%
%   To create a GaussNoise object, one of two constructors may be used.
%   The default constructor, GaussNoise(), will set all properties to their
%   default values.  The alternative constructor adopts MATLAB's
%   property/value string pair convention to allow the user to specify any
%   subset of properties, with remaining properties set to their default
%   values.  For example, GaussNoise('prior_var', 1e-3, 'learn_prior_var', 
%   'column') will construct a GaussNoise object that specifies an initial
%   variance PSI(m,t) = 1e-3 for all m,t, and indicates that a unique
%   variance is to be learned using the EM algorithm for each column of W.
%
%   Note also that if the user wishes to have EMTurboGAMP attempt to 
%   intelligently initialize PSI based on the matrix of measurements, Y, 
%   and the sensing matrix A, then GaussNoise('init_params', 'true') is the
%   appropriate constructor.
%
% PROPERTIES (State variables)
%   prior_var           The prior variance, PSI.  [Default: 1e-2]
%   learn_prior_var     Specifies options for learning the noise variance
%                       using an EM learning algorithm.  Acceptable
%                       settings are 'scalar', 'row', 'column', or 'false'.
%                       [Default: 'scalar']
%   init_params     	Initialize PSI automatically using Y and A 
%                       ('true'), or not ('false').  [Default: 'false']
%   data                Noise is real-valued ('real') or complex-valued
%                       ('complex')? [Default: 'real']
%   version             Execute sum-product GAMP for MMSE estimation
%                       ('mmse') or max-sum GAMP for MAP estimation ('map')
%                       [Default: 'mmse']
%
% METHODS (Subroutines/functions)
%   GaussNoise()
%       - Default constructor.  Assigns all properties to their default 
%         values.
%   GaussNoise('ParameterName1', Value1, 'ParameterName2', Value2)
%       - Custom constructor.  Can be used to set any subset of the
%         parameters, with remaining parameters set to their defaults
%   set(obj, 'ParameterName', Value)
%       - Method that can be used to set the value of a parameter once the
%         object of class GaussNoise, obj, has been constructed
%   print(obj)
%       - Print the current value of each property to the command window
%   get_type(obj)
%       - Returns the character string 'AWGN' when obj is a GaussNoise
%         object
%   GaussNoiseCopyObj = copy(obj)
%       - Creates an independent copy of a GaussNoise object, obj
%   EstimOut = UpdatePriors(obj, GAMPState, Y, A)
%     	- Given the final state of the message passing variables that were
%         output from GAMP after its most recent execution, produce a new 
%         object of the EstimOut base class that will be used to specify 
%         the AWGN noise "prior" on the next iteration of GAMP. TBobj is an 
%         object of the TurboOpt class, GAMPState is an object of the 
%         GAMPState class, and Y is an M-by-T matrix of measurements.  If 
%         TBobj.commonA is false, then this method returns a 1-by-T cell
%         array of EstimOut objects.  [Hidden method]
%   EstimOut = InitPriors(TBobj, Y)
%    	- Provides an initial EstimOut object for use by GAMP the first
%         time. If the user wishes to initialize parameters from the
%         data, Y, then argument Y must be provided.  TBobj is a TurboOpt
%         object.  If TBobj.commonA is false, then this method returns a 
%         1-by-T cell array of EstimOut objects. [Hidden method]
%   NOISE = genRand(TBobj, GenParams)
%       - Produce a realization of AWGN noise matrix, given a TurboOpt
%         object (TBobj) and GenParams object as inputs [Hidden method]
%

%
% Coded by: Justin Ziniel, The Ohio State Univ.
% E-mail: zinielj@ece.osu.edu
% Last change: 09/18/12
% Change summary: 
%       - Created (10/14/11; JAZ)
%       - Added implementation of UpdatePriors, which is called by the main
%         EMturboGAMP function to produce new EstimOut object (12/13/11;
%         JAZ)
%       - Added genRand method (01/03/12; JAZ)
%       - Modified description for a time-varying A(t) model (02/07/12;
%         JAZ)
%       - Added support for complex-valued noise (05/22/12; JAZ)
%       - Added copy method (07/10/12; JAZ)
%       - Added support for max-sum GAMP operation (09/18/12; JAZ)
% Version 0.2
%

classdef GaussNoise < Observation
    
    properties
        prior_var = 1e-2;               % Prior noise variance
        learn_prior_var = 'scalar';     % Learn a single noise variance
        init_params = 'false';          % Don't initialize params from data
        data = 'real';                  % Real- or complex-valued noise?
        version = 'mmse';               % Sum-product ('mmse') or max-sum ('map') GAMP
    end % properties
    
    properties (Access = private, Hidden)
        InitIter = true;        % Flag to prevent EM learning on 1st 
                                % smoothing iteration
    end
    
    properties (Constant, Hidden)
        type = 'AWGN';   % Noise distribution identifier string
    end
   
    methods 
        
        % *****************************************************************
        %                      CONSTRUCTOR METHOD
        % *****************************************************************
        
        function obj = GaussNoise(varargin)
            if nargin == 1 || mod(nargin, 2) ~= 0
                error('Improper constructor call')
            else
                for i = 1 : 2 : nargin - 1
                    obj.set(lower(varargin{i}), varargin{i+1});
                end
            end
        end
        
        
        % *****************************************************************
        %                         SET METHODS
        % *****************************************************************
        
        % Set method for prior_var (PSI)
        function obj = set.prior_var(obj, PSI)
           if min(min(PSI)) < 0
              error('Noise variances must be non-negative')
           else
              obj.prior_var = PSI;
           end
        end
        
        % Set method for learn_prior_var
        function obj = set.learn_prior_var(obj, string)
            if ~check_input(obj, string)
                error('Invalid option: learn_prior_var')
            end
            obj.learn_prior_var = lower(string);
        end
        
        % Set method for init_params
        function obj = set.init_params(obj, string)
            if sum(strcmpi(string, {'true', 'false'})) == 0
                error('Invalid option: init_params')
            end
            obj.init_params = lower(string);
        end
        
        % Set method for data
        function obj = set.data(obj, string)
            if sum(strcmpi(string, {'real', 'complex'})) == 0
                error('Invalid option: data')
            end
            obj.data = lower(string);
        end
        
        % Set method for version
        function obj = set.version(obj, version)
            if any(strcmpi(version, {'mmse', 'map'}))
                obj.version = lower(version);
            else
                error('Invalid option: version')
            end
        end
        
        
        % *****************************************************************
        %                         PRINT METHOD
        % *****************************************************************
        
        % Print the current configuration to the command window
        function print(obj)
            fprintf('NOISE PRIOR: AWGN\n')
            fprintf('             prior_var: %s\n', ...
                form(obj, obj.prior_var))
            fprintf('       learn_prior_var: %s\n', obj.learn_prior_var)
            fprintf('           init_params: %s\n', obj.init_params)
            fprintf('                  data: %s\n', obj.data)
            fprintf('          GAMP version: %s\n', obj.version)
        end
        
        
        % *****************************************************************
        %                         COPY METHOD
        % *****************************************************************
        
        % Creates an independent copy of a GaussNoise object
        function GaussNoiseCopyObj = copy(obj)
            GaussNoiseCopyObj = GaussNoise('prior_var', obj.prior_var, ...
                'learn_prior_var', obj.learn_prior_var, 'init_params', ...
                obj.init_params, 'data', obj.data, 'version', obj.version);
        end
        
        % *****************************************************************
        %                        ACCESSORY METHOD
        % *****************************************************************
        
        % This function allows one to query which type of noise family
        % this class is by returning a character string type identifier
        function type_id = get_type(obj)
            type_id = obj.type;
        end
        
    end % methods
    
    methods (Hidden)
        
        % *****************************************************************
        %      	         UPDATE GAMP NOISE "PRIOR" METHOD
        % *****************************************************************
        
      	% Update the EstimOut object for an AWGN noise prior
        function EstimOut = UpdatePriors(obj, TBobj, GAMPState, Y, ...
                ~, ~)
            
            [M, T] = size(Y);
            
            % Unpack the GAMPState object
            [~, ~, ~, ~, ~, ~, Zhat, Zvar] = GAMPState.getState();
            
            % Perform EM update of AWGN variance, if user wishes.  Note: We
            % ignore the addition of sum(Zvar)/M/T to these updates, as it
            % appears to cause convergence issues
            if strcmpi('scalar', obj.learn_prior_var)
                % Learn a scalar noise variance
                var_upd = sum(sum(abs(Y - Zhat).^2 + Zvar)) / M / T;
            elseif strcmpi('row', obj.learn_prior_var)
                % Learn a unique noise variance for each row
                var_upd = sum(abs(Y - Zhat).^2 + Zvar, 2) / T;
            elseif strcmpi('column', obj.learn_prior_var)
                % Learn a unique noise variance for each column
                var_upd = sum(abs(Y - Zhat).^2 + Zvar, 1) / M;
            else
                % No refinement
                var_upd = obj.prior_var;
            end
            
            % If this is the first smoothing iteration, avoid attempting to
            % refine the noise variance
            if obj.InitIter
                var_upd = obj.prior_var;    % Use existing noise variance
                obj.InitIter = false;       % Clear the flag
            end
                
            
            % Place an M-by-T matrix of updated variances back into object
            % structure
            obj.prior_var = TBobj.resize(var_upd, M, T);
            
            % Identify whether to use sum-product or max-sum GAMP EstimOut
            % objects
            maxSumVal = strcmpi(obj.version, 'map');
            
            % Construct an AwgnEstimOut object
            switch TBobj.commonA
                case true
                    % There is one common A matrix for all timesteps, thus
                    % we can use a matrix-valued EstimOut object and run
                    % matrix GAMP
                    switch obj.data
                        case 'real'
                            EstimOut = AwgnEstimOut(Y, obj.prior_var, ...
                                maxSumVal);
                        case 'complex'
                            EstimOut = CAwgnEstimOut(Y, obj.prior_var, ...
                                maxSumVal);
                    end
                case false
                    % Each timestep has a different matrix A(t), thus we
                    % must run vector GAMP with separate EstimOut objects
                    % for each timestep.  Return a cell array of such
                    % objects
                    EstimOut = cell(1,T);
                    for t = 1:T
                        switch obj.data
                            case 'real'
                                EstimOut{t} = AwgnEstimOut(Y(:,t), ...
                                    obj.prior_var(:,t), maxSumVal);
                            case 'complex'
                                EstimOut{t} = CAwgnEstimOut(Y(:,t), ...
                                    obj.prior_var(:,t), maxSumVal);
                        end
                    end
            end
        end
        
        
        % *****************************************************************
        %          	   INITIALIZE GAMP NOISE "PRIOR" METHOD
        % *****************************************************************
        
        % Initialize EstimOut object for AWGN noise prior
        function EstimOut = InitPriors(obj, TBobj, Y)
            
            % First get the problem size
            [M, T] = size(Y);
            
            % Initialize differently based on user preference
            switch TBobj.Observation.init_params
                case 'false'
                    % Resize user-provided variance to make M-by-T
                    TBobj.Observation.prior_var = ...
                        TBobj.resize(TBobj.Observation.prior_var, M, T);
                    
                case 'true'
                    % Make a crude guess of the noise variance, psi
                    if strcmp(obj.learn_prior_var, 'column')
                        psi = sum(abs(Y).^2, 1) / M / 101;
                    else
                        psi = norm(Y, 'fro')^2 / (M*T*101);
                    end
                    
                    TBobj.Observation.prior_var = TBobj.resize(psi, M, T);
            end
            
            % Identify whether to use sum-product or max-sum GAMP EstimOut
            % objects
            maxSumVal = strcmpi(obj.version, 'map');
            
            % Form EstimOut object
            switch TBobj.commonA
                case true
                    % There is one common A matrix for all timesteps, thus
                    % we can use a matrix-valued EstimOut object and run
                    % matrix GAMP
                    switch obj.data
                        case 'real'
                            EstimOut = AwgnEstimOut(Y, obj.prior_var, ...
                                maxSumVal);
                        case 'complex'
                            EstimOut = CAwgnEstimOut(Y, obj.prior_var, ...
                                maxSumVal);
                    end
                case false
                    % Each timestep has a different matrix A(t), thus we
                    % must run vector GAMP with separate EstimOut objects
                    % for each timestep.  Return a cell array of such
                    % objects
                    EstimOut = cell(1,T);
                    for t = 1:T
                        switch obj.data
                            case 'real'
                                EstimOut{t} = AwgnEstimOut(Y(:,t), ...
                                    obj.prior_var(:,t), maxSumVal);
                            case 'complex'
                                EstimOut{t} = CAwgnEstimOut(Y(:,t), ...
                                    obj.prior_var(:,t), maxSumVal);
                        end
                    end
            end
        end
        
        % *****************************************************************
        %                   GENERATE NOISE REALIZATION
        % *****************************************************************
        
        function NOISE = genRand(obj, TBobj, GenParams)
            % Get size information
            M = GenParams.M;
            T = GenParams.T;
            
            noisevar = TBobj.resize(obj.prior_var, M, T);
            
            % Generate a realization
            switch GenParams.type
                case 'real'
                    NOISE = sqrt(noisevar) .* randn(M,T);
                case 'complex'
                    NOISE = sqrt(noisevar/2) .* randn(M,T) + ...
                        1j*sqrt(noisevar/2) .* randn(M,T);
            end
        end
    end
    
    methods (Access = private)
        
        % *****************************************************************
        %                        HELPER METHODS
        % *****************************************************************
        
        % This method makes sure that the inputs to learn_sparsity_rate,
        % etc. are valid options
        function flag = check_input(obj, string)
            flag = 1;
            if ~isa(string, 'char')
                flag = 0;
            elseif sum(strcmpi(string, {'scalar', 'row', 'column', ...
                    'false'})) == 0
                flag = 0;
            end
        end
        
        % This method is called by the print method in order to format how
        % properties are printed.  Properties that are scalars will have
        % their values printed to the command window, while arrays will
        % have their dimensions, minimal, and maximal values printed to the
        % screen
        function string = form(obj, prop)
            if numel(prop) == 1
                string = num2str(prop);
            else
                string = sprintf('%d-by-%d array (Min: %g, Max: %g)', ...
                    size(prop, 1), size(prop, 2), min(min(prop)), ...
                    max(max(prop)));
            end
        end
    end % methods
    
   
end % classdef
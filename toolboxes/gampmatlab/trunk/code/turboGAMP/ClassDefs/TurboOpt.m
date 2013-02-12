% CLASS: TurboOpt
% 
% HIERARCHY (Enumeration of the various super- and subclasses)
%   Superclasses: hgsetget (MATLAB handle class)
%   Subclasses: N/A
% 
% TYPE (Abstract or Concrete)
%   Concrete
%
% DESCRIPTION (High-level overview of the class)
%   The TurboOpt class contains properties that collectively define the
%   setup of a particular compressive sensing (CS) recovery problem.
%   Specifically, the properties (a.k.a. variables) within the TurboOpt
%   class describe the probabilistic model underlying the sparse linear
%   inverse problem in which a length-N unknown sparse signal vector x(t)
%   is transformed by a linear operator A(t) into a length-M "transform
%   vector" z(t), i.e.,
%                    z(t) = A(t)*x(t),    t = 1, ..., T   
%   which is observed through a separable noisy channel, yielding a
%   length-M vector of real- or complex-valued observations y(t).  The
%   relationship between the transform vector, z(t), and the observation
%   vector, y(t), is described by the separable probability distribution
%       p(y(t)|z(t)) = p(y(1,t)|z(1,t)) * ... * p(y(M,t)|z(M,t)).
%   
%   Note that there are several special cases of this signal model.  The
%   first is the standard time-varying CS measurement model:
%               	y(t) = A(t)*x(t) + w(t),  t = 1, ..., T
%   where y(t) is a length-M vector of observations, A(t) is an M-by-N
%   sensing matrix, x(t) is a length-N unknown signal vector, and w(t) is 
%   a length-M vector of corrupting Observation.  Another special case is  
%   the MMV signal model
%                             Y = A*X + W.
%   Lastly, the single measurement vector (SMV) CS recovery problem is
%   a special case as well, (i.e., T = 1).
%
%   The two most important properties of the TurboOpt class are the Signal
%   and Observation properties.  The Signal property is an object that is a
%   subclass of the abstract Signal class (see Signal.m), e.g., a BernGauss
%   object (see BernGauss.m).  It contains the parameters that are needed
%   to define the marginal prior distributions of the elements of X.
%   Likewise, the Observation property is an object that is a subclass of 
%   the abstract Observation class (see Observation.m), e.g., a GaussNoise 
%   object (see GaussNoise.m).  The Observation property contains the 
%   parameters needed to define the channel Observation distribution, 
%   p(y(t)|z(t)) (e.g., distribution of W).
%
%   There are multiple valid ways to construct a TurboOpt object, which we
%   will call TBobj.  If the user wishes to set up a CS problem in which
%   the elements of X are i.i.d. Bernoulli-Gaussian, the Observation is 
%   AWGN, and there is no model of structured sparsity, the default 
%   constructor TurboOpt() can be used, (e.g., TBobj = TurboOpt();).  Once 
%   TBobj has been constructed, it is possible to modify specific signal 
%   model parameters in the same way one manipulates structures in MATLAB.  
%   For instance, to change the prior sparsity rate of the Bernoulli-Gaussian
%   prior on X, one simply executes TBobj.Signal.sparsity_rate = LAMBDA, 
%   where LAMBDA is a valid assignment, as defined by the BernGauss class.
%   Alternatively, one can construct TBobj by providing the constructor
%   with valid objects derived from the corresponding base classes.  For
%   instance, suppose BGobj is an object of the BernGauss class that has
%   already been constructed and initialized to yield a desired Bernoulli-
%   Gaussian prior, and suppose AWGNobj is a similarly constructed and
%   initialized GaussNoise object.  Since BernGauss and GaussNoise are
%   subclasses of the abstract Signal and Observation classes, respectively,
%   BGobj and AWGNobj constitute valid assignments to the Signal and 
%   Observation properties of a TurboOpt object.  Therefore, the 
%   constructor TBobj = TurboOpt('Signal', BGobj, 'Observation', AWGNobj) 
%   can be used to create a TurboOpt object with the desired signal model.  
%   More generally, any subset of TurboOpt's properties can be set in the 
%   constructor by following MATLAB's property/value string pair 
%   convention, with remaining properties set to their defaults.
%
% PROPERTIES (State variables)
%   Signal          An object that is a concrete subclass of the abstract
%                   Signal class, used to define the marginal priors of the
%                   elements of X.  [Default: A default BernGauss object]
%   Observation    	An object that is a concrete subclass of the abstract
%                   Observation class, used to define the channel model,
%                   p(Y | Z).  [Default: A default GaussNoise  object]
%   SupportStruct   An object of the SupportStruct class, or of an
%                   inheriting sub-class (e.g., MarkovChain1).  This
%                   property is used to specify the form of structure that
%                   is found in the support pattern of the signal X.
%                   [Default: A default SupportStruct object (support is
%                   independent)]
%   AmplitudeStruct An object of the AmplitudeStruct class, or of an
%                   inheriting sub-class (e.g., GaussMarkov).  This
%                   property is used to specify the form of structure that
%                   is found in the amplitudes of the non-zero elements of
%                   the signal X. [Default: A default AmplitudeStruct
%                   object (amplitudes are independently distributed
%                   according to the distribution of active elements
%                   specified in the Signal class)]
%   RunOptions      An object of the RunOptions class, which contains
%                   runtime parameters governing the execution of
%                   EMturboGAMP.  [Default: A default RunOptions object]
%
% METHODS (Subroutines/functions)
%   TurboOpt()
%       - Default constructor.  Creates a default signal model in which the
%         elements of X are i.i.d. Bernoulli-Gaussian random variables,
%         with a distribution specified by the default values of a
%         BernGauss object, and in which the Observation W is modelled as 
%         AWGN with statistics given by the defaults of a GaussNoise 
%         object.  No assumption is made with regards to structured 
%         sparsity.
%   TurboOpt('PropertyName1', Value1, 'PropertyName2', Value2, ...)
%       - Custom constructor.  Can be used to set any subset of the
%         properties, with remaining properties set to their defaults
%   [EstimIn, EstimOut, S_POST] = UpdatePriors(obj, GAMPState, Y, ...
%       EstimInOld, EstimOutOld, A)
%       - Produces objects derived from the EstimIn and EstimOut base
%         classes used by GAMP, based on the outputs from the previous
%         execution of GAMP (contained in GAMPState object), and the M-by-T
%         matrix of measurements, Y, along with an estimated posterior for
%         the support variables in S_POST
%   [EstimIn, EstimOut] = InitPriors(obj, Y, A)
%       - Initializes objects derived from the EstimIn and EstimOut base
%         classes, which can then be used for an initial iteration of the
%         GAMP algorithm.
%   set(obj, 'PropertyName', Value)
%       - Method that can be used to set the value of a parameter once the
%         object of class TurboOpt, obj, has been constructed
%   print(obj)
%       - Print the current value of each property of obj to the command 
%         window
%   TurboCopyObj = copy(obj)
%       - Since TurboOpt is a handle class (inheriting from hgsetget), in
%         order to create an independent copy of a TurboOpt object (i.e.,
%         one that does not point to the same underlying object), call this
%         copy method to create an independent copy of the TurboOpt object,
%         obj, in TurboCopyObj
%   delete(obj)
%       - Destroy the TurboOpt object, obj.
%       
%

%
% Coded by: Justin Ziniel, The Ohio State Univ.
% E-mail: zinielj@ece.osu.edu
% Last change: 07/10/12
% Change summary: 
%       - Created (10/16/11; JAZ)
%       - Added methods UpdatePriors and InitPriors, and AmplitudeStruct
%         property (12/15/11; JAZ)
%       - Added RunOptions property, generalized the signal model and
%         included support for time-varying linear transform operators, 
%         A(t) (02/07/12; JAZ)
%       - Added S_POST as a variable returned by UpdatePriors method
%         (04/12/12; JAZ)
%       - Renamed Observation class references to Observation class 
%         (05/22/12; JAZ)
%       - Created a copy method, useful for obtaining an independent clone
%         of a TurboOpt object (07/10/12; JAZ)
% Version 0.2
%

classdef TurboOpt < hgsetget
    
    properties 
        Signal;
        Observation;
        SupportStruct;
        AmplitudeStruct;
        RunOptions;
    end % properties
    
    properties (Hidden)
        commonA = true;  % Logical that denotes whether A(t) = A for all t
    end
   
    methods
                
        % *****************************************************************
        %                      CONSTRUCTOR METHOD
        % *****************************************************************
        
        % Multi-argument constructor
        function obj = TurboOpt(varargin)
            if nargin == 1 || mod(nargin, 2) ~= 0
                error('Improper constructor call')
            else
                % First make a default TurboOpt object, then fill in
                % user-provided values
                obj.Signal = BernGauss();   % Default signal prior is B.G.
                obj.Observation = GaussNoise();   % Default observation prior is AWGN
                obj.SupportStruct = NoSupportStruct();	% Default support is independent
                obj.AmplitudeStruct = NoAmplitudeStruct();    % No amplitude structure
                obj.RunOptions = RunOptions();  % Default runtime parameters
                for i = 1 : 2 : nargin - 1
                    obj.set(varargin{i}, varargin{i+1});
                end
            end
        end
        
        
        % *****************************************************************
        %                    GAMP INTERFACE METHODS
        % *****************************************************************
        
        % UPDATEPRIORS - Produce new EstimIn and EstimOut objects for GAMP
        % to use in its recoveries, based on messages received from the
        % last round of recoveries
        %
        % INPUTS:
        % obj           Object of the TurboOpt type
        % GAMPState     GAMPState object, containing final state of GAMP
        % Y             M-by-T matrix of measurements
        % EstimInOld    The most recent EstimIn object given to GAMP
        % EstimOutOld   The most recent EstimOut object given to GAMP
        % A             Linear operator, or 1-by-T cell array of linear
        %               operators
        %
        % OUTPUTS:
        % EstimIn       Concrete subclass of the EstimIn class
        % EstimOut      Concrete subclass of the EstimOut class
        function [EstimIn, EstimOut, S_POST] = UpdatePriors(obj, ...
                GAMPState, Y, EstimInOld, EstimOutOld, A)
            % First do some basic input checking
            if nargin < 3
                error('Insufficient number of input arguments')
            end
            if ~isa(obj, 'TurboOpt')
                error('UpdatePriors requires a TurboOpt object')
            end
            if ~isa(GAMPState, 'GAMPState')
                error('UpdatePriors requires a GAMPState object')
            end
            
            % First use GAMP's final state to generate a new EstimIn
            % object
            SigObj = obj.Signal;
            [EstimIn, S_POST] = ...
                SigObj.UpdatePriors(obj, GAMPState, EstimInOld);
            
            % Now use GAMP's final state to generate a new EstimOut object
            ObservationObj = obj.Observation;
            EstimOut = ObservationObj.UpdatePriors(obj, GAMPState, Y, ...
                EstimOutOld, A);
        end
        
        
        % INITPRIORS - Produce initial EstimIn and EstimOut objects for 
        % GAMP to use for its first iteration
        %
        % INPUTS:
        % obj           Object of the TurboOpt type
        % Y             M-by-T matrix of measurements
        % A             M-by-N sensing matrix, or an object derived from
        %               GAMP's LinTrans class
        %
        % OUTPUTS:
        % EstimIn       Concrete subclass of the EstimIn class
        % EstimOut      Concrete subclass of the EstimOut class
        function [EstimIn, EstimOut] = InitPriors(obj, Y, A)
            % First do some basic input checking
            if nargin < 3
                error('Insufficient number of input arguments')
            end
            if ~isa(obj, 'TurboOpt')
                error('InitPriors requires a TurboOpt object')
            end
            if ~isa(Y, 'double')
                error(sprintf('%s must be an array of measurements\n', ...
                    inputname(2)))
            end
            if ~isa(A, 'double') && ~isa(A, 'LinTrans') && ~isa(A, 'cell')
                error(['A must be a 1-by-T cell array, or a single ' ...
                    'explicit matrix, or a subclass of GAMP''s LinTrans ' ...
                    'class'])
            end
            
            % If A is a cell array, we will not be using matrix-valued GAMP
            if isa(A, 'cell')
                obj.commonA = false;
            end
            
            % Produce an initial EstimIn object
%             SigObj = obj.Signal;
%             EstimIn = SigObj.InitPriors(obj, Y, A);
            EstimIn = obj.Signal.InitPriors(obj, Y, A);
            
            % Produce an initial EstimOut object
%             ObservationObj = obj.Observation;
%             EstimOut = ObservationObj.InitPriors(obj, Y);
            EstimOut = obj.Observation.InitPriors(obj, Y);

            % Now switch based on the class of amplitude structure in order
            % to initialize any relevant amplitude structure parameters
            switch obj.AmplitudeStruct.get_type()
                case 'GM'   % Gauss-Markov amplitude structure
                    AmplObj = obj.AmplitudeStruct;
                    AmplObj.InitPriors(obj, Y, A);
                otherwise
                    % Nothing to call here
            end
        end
        
        
        % *****************************************************************
        %                       SET METHODS
        % *****************************************************************
        
        % Set method for Signal property
        function obj = set.Signal(obj, SigObj)
            % First verify that the input SigObj is a valid subclass of the
            % Signal superclass
            if ~isa(SigObj, 'Signal')
                error('TurboOpt.Signal must be a valid Signal subclass')
            else
                obj.Signal = SigObj;
            end
        end
        
        % Set method for Observation property
        function obj = set.Observation(obj, ObservationObj)
            % First verify that the input NoiseObj is a valid subclass of 
            % the Observation superclass
            if ~isa(ObservationObj, 'Observation')
                error('TurboOpt.Observation must be a valid Observation subclass')
            else
                obj.Observation = ObservationObj;
            end
        end
        
        % Set method for SupportStruct property
        function obj = set.SupportStruct(obj, SuppObj)
            % First verify that the input SuppObj is a valid object of the
            % SupportStruct class, or of an inheriting subclass
            if ~isa(SuppObj, 'SupportStruct')
                error(['TurboOpt.SupportStruct must be a valid ' ...
                'SupportStruct or inheriting subclass object'])
            else
                obj.SupportStruct = SuppObj;
            end
        end
        
        % Set method for AmplitudeStruct property
        function obj = set.AmplitudeStruct(obj, AmplObj)
            % First verify that the input SuppObj is a valid object of the
            % SupportStruct class, or of an inheriting subclass
            if ~isa(AmplObj, 'AmplitudeStruct')
                error(['TurboOpt.AmplitudeStruct must be a valid ' ...
                'AmplitudeStruct or inheriting subclass object'])
            else
                obj.AmplitudeStruct = AmplObj;
            end
        end
        
        % Set method for AmplitudeStruct property
        function obj = set.RunOptions(obj, RunObj)
            % First verify that the input SuppObj is a valid object of the
            % SupportStruct class, or of an inheriting subclass
            if ~isa(RunObj, 'RunOptions')
                error(['TurboOpt.RunOptions must be a valid ' ...
                'RunOptions object'])
            else
                obj.RunOptions = RunObj;
            end
        end
        
                
        % *****************************************************************
        %                           PRINT METHOD
        % *****************************************************************
        
        % Print the value of each property to the command window, using the
        % print methods of each class
        function print(obj)
            obj.Signal.print();
            obj.Observation.print();
            obj.SupportStruct.print();
            obj.AmplitudeStruct.print();
            obj.RunOptions.print();
        end
        
        
        % *****************************************************************
        %                           COPY METHOD
        % *****************************************************************
        
        % Create an independent clone of a TurboOpt object
        function TurboCopyObj = copy(obj)
            TurboCopyObj = TurboOpt('Signal', obj.Signal.copy(), ...
                'Observation', obj.Observation.copy(), 'SupportStruct', ...
                obj.SupportStruct.copy(), 'AmplitudeStruct', ...
                obj.AmplitudeStruct.copy(), 'RunOptions', ...
                obj.RunOptions.copy());
        end
        
        
        % *****************************************************************
        %                          DELETE METHOD
        % *****************************************************************
        
        % Destroy all the handles to each of the properties
        function delete(obj)
            obj.Signal.delete();
            obj.Observation.delete();
            obj.SupportStruct.delete();
            obj.AmplitudeStruct.delete();
            obj.RunOptions.delete();
        end
    end % methods
    
            
    methods (Hidden)
        % *****************************************************************
        %                         HELPER METHODS
        % *****************************************************************
        % RESIZE: A function that will resize input IN to make it have
        % dimension J-by-K-by-L, by replicating it if need be
        function OUT = resize(obj, IN, J, K, L)
            if nargin < 5
                % Assume third dimension is singleton
                L = 1;
            end
            if size(IN, 1) == 1 && size(IN, 2) == 1 && size(IN, 3) == L
                OUT = repmat(IN, [J, K, 1]);
            elseif size(IN, 1) == J && size(IN, 2) == 1 && size(IN, 3) == L
                OUT = repmat(IN, [1, K, 1]);
            elseif size(IN, 1) == 1 && size(IN, 2) == K && size(IN, 3) == L
                OUT = repmat(IN, [J, 1, 1]);
            elseif size(IN, 1) == 1 && size(IN, 2) == K && size(IN, 3) == 1
                OUT = repmat(IN, [J, 1, L]);
            elseif size(IN, 1) == J && size(IN, 2) == K && size(IN, 3) == L
                OUT = IN;
            elseif size(IN, 1) == 1 && size(IN, 2) == 1 && size(IN, 3) == 1
                OUT = repmat(IN, [J, K, L]);
            else
                error('Incorrect size of parameter')
            end
        end
    end
   
end % classdef
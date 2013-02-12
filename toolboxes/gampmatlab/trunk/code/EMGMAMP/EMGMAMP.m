%EMGMAMP  Expectation-Maximization Gaussian-Mixture AMP
% From the (possibly complex-valued) M-by-T matrix of noisy linear 
% observations Y = AX + W, EMGMAMP returns approximately MMSE estimates 
% of the N-by-T signal matrix X, where the focus is the sub-Nyquist case 
% M < N, as in sparse reconstruction and compressive sensing.  In EMGMAMP,
% it is assumed that the signal coefficients X(n,t) are apriori independent 
% and Bernoulli-Gaussian-mixture distributed according to the pdf:
%
%   p(X(n,t)) = (1-lambda(t))*delta(X(n,t))
%      + sum_k^L lambda(t)*omega(t,k) (C)Normal(X(n,t);theta(t,k),phi(t,k))
%
% where n = 1:N indexes the row of the signal matrix X,  t = 1:T 
% indexes the column (timestep) of all matrices, and delta() denotes the
% Dirac delta pdf.  Above, lambda(t), omega(t,:), theta(t,:), and phi(t,:) 
% are the sparsity rate, mixture weights, active-coefficient means, and 
% active-coefficient variances, respectively, of the coefficients within 
% the t-th column of X.  In EMGMAMP, the noise samples W(m,t) are assumed 
% apriori independent and zero-mean Gaussian distributed: 
%
%    p(W(m,t)) = (C)Normal(W(m,t);0,psi(t))
%
% where m = 1:M indexes the rows of the matrices W and Y, and psi(t) 
% is the variance within the t-th column of W.
%
% The estimate of X is calculated using the Generalized Approximate Message 
% Passing (GAMP) algorithm, as described in "Generalized Approximate Message 
% Passing for Estimation with Random Linear Mixing" by Sundeep Rangan,
% with certain adaptive-stepsize enhancements.
%
% The parameters lambda, omega, theta, phi, and psi are automatically 
% learned via the Expectation Maximization (EM) algorithm, using X and W as 
% the "hidden" data. 
%
% The EMGMAMP algorithm is described in detail in "Expectation-Maximization
% Gaussian-Mixture Approximate Message Passing" by Jeremy Vila and 
% Philip Schniter.
%
%Syntax:
% [Xhat, Xvar, param] = EMGMAMP(Y, A) % suggested default 
% [Xhat, Xvar, param] = EMGMAMP(Y, A, GAMPopt, EMopt) % allows customization 
%
%Inputs:
% Y - matrix of observed measurements
%
% A - known mixing matrix (in either explicit or GAMP-object format)
%
% GAMPopt - structure containing various GAMP option fields [optional].
%           User MUST call the constructor first (i.e., GAMPopt =
%           GampOpt())
%   .nit                number of iterations
%   .pvarMin            minimum variance of each element of p
%   .XvarMin            minimum variance of each element of x
%   .step               step size 
%   .stepMin            minimum step size 
%   .stepMax            maximum step size
%   .stepIncr           Multiplicative step size increase, when successful
%   .stepDecr           Multiplicative step size decrease, when unsuccessful
%   .pvarStep           Logical flag to include a step size in the pvar/zvar calculation.
%   .varNorm            %Option to "normalize" variances for computation.
%   .adaptStep          adaptive step size [important if A column-correlated]
%   .stepWindow         step size check window size
%   .bbStep             Barzilai Borwein step size [if A column-correlated]
%   .verbose            Print results in each iteration
%   .tol                Convergence tolerance
%   .stepTol            minimum allowed step size
%   .Avar               variance in A entries (may be scalar or matrix)
%   .xhat0              the initialization of x for GAMP [default handled
%                       internally]
%   .xvar0              the intialization of variance of each x_n for GAMP
%                       [default handled internally]
%   .shat0              the initialization of s for GAMP [default handled
%                       internally]
%
%   Note: EMGMAMP(Y, A) or EMGMAMP(Y, A, []) sets GAMP parameters at defaults
%
% EMopt - structure containing various EM option fields [optional]
%   .heavy_tailed       Set to true for heavy-tailed (compressible or
%                       non-compressible) signals.  [default=false]
%   .maxEMiter          number of EM iterations [default=20]
%   .EMtol              EMGMAMP terminates when norm change between 
%                       sparsity rate is less than this tolerance. Set to
%                       -1 to turn off and run max iters [default=1e-5]
%   .learn_lambda       Set to true to learn lambda, set to false
%                       to never update lambda (stays at initialization) 
%                       [default=true]
%   .learn_weights      Set to true to learn the weights omega of the 
%                       Gaussian Mixture [default=true]
%   .learn_mean         Set to true to learn active mean (theta), set to false 
%                       to never update theta (stays at initialization) 
%                       [default=true]
%   .learn_var          Set to true to learn active variance (phi), set to 
%                       false to never update active variance (stays at  
%                       initialization) [default=true]
%   .learn_noisevar     Set to true to learn noise variance (psi), set to
%                       false to never update noise variance (stays at  
%                       initialization) [default=true]
%   .sig_dim            Set to 'col' to allow different EM-learned parameters
%                       for each column in signal matrix X.  Set to 'joint' 
%                       to force common EM-learned parameters for all columns.
%                       [default='col']
%
%   Note: EMGMAMP(Y, A) or EMGMAMP(Y, A, [], []) sets EM params at defaults
%
%   WARNING! Set the following initializations only if confident they
%   accurately represent the model.  Otherwise, let the defaults preside.
%   Also, ensure that all GM-component vectors are of compatible lengths.  
%
%   .L                  Number of mixture components [default=3, max=20]
%   .lambda             Initialize overall sparsity rate lambda. Must be 
%                       either scalar or 1-by-T vector.  [default is 
%                       based on the LASSO phase transition curve]
%   .active_weights     Initialize GM component weights (omega).
%                       Must be L-by-1 or L-by-T vector.
%                       [Defaults based on best L-term uniform-pdf fit]
%   .active_mean        Initialize mean of active components (theta)
%                       Must be L-by-1 or L-by-T vector.
%                       [Defaults based on best L-term uniform-pdf fit]
%   .active_var         Initialize variance of active components (phi)
%                       Must be L-by-1 or L-by-T vector.
%                       [Defaults based on best L-term uniform-pdf fit]
%   .noise_var          Initialize the noise variance.  Must be either 
%                       scalar or 1-by-T vector [default is based on an 
%                       assumed SNR = 20 dB]
%   
%   Note: Setting EMopt.heavy_tailed=true overrides the above parameter 
%   initializations, and sets EMopt.learn_mean=false.
%
%Outputs:
% Xhat - GAMP-estimated posterior means
% Xvar - GAMP-estimated posterior variances
% param - structure containing EM-estimated prior parameters
%   .lambda             Sparsity rate (lambda) [1-by-T for 'col' mode]
%   .active_weights     weights of GM (omega) [L-by-T for 'col' mode]
%   .active_mean        Means of active components (theta) [L-by-T for 'col' mode]
%   .active_var         Variances of active components (phi) [L-by-T for 'col' mode]
%   .noise_var          Variance of noise (psi) [1-by-T for 'col' mode]
% estHist - GAMP output structure containing per-iteration GAMP data
%   .pass               Per iteration logical value indicating whether 
%                       adaptive stepsize selection was used.
%   .step               Per-iteration step sizes.
%   .val                Per-iteration value of the GAMP cost function.
% 
%
% Coded by: Jeremy Vila, The Ohio State Univ.
% E-mail: vilaj@ece.osu.edu
% Last change: 9/04/12
% Change summary: 
%   v 1.0 (JV)- First release
%   v 1.1 (JV)- Fixed some minor bugs
%   v 2.0 (JV)- Accounts for MMV model.  Also made changes due to changes
%               in GAMP.
%   v 2.1 (JV)- Changed the way GAMPopt is initialized.  Use must call
%               constructor first.
%   v 2.2 (JV)- Updated noise variance learning in accordance to new GAMP
%               outputs
%   v 3.0 (JV)- Added option to output GAMP Histories.
%
% Version 3.0
%
function [Xhat, Xvar, param, estHist] = EMGMAMP(Y, A, GAMPopt, EMopt)

%Ensure that matrix is GAMP object
if ~isobject(A)
    A = MatrixLinTrans(A); 
end

%Find problem dimensions
[M,N] = A.size();
T = size(Y, 2);
if size(Y, 1) ~= M
    error('Dimension mismatch betweeen Y and A')
end

%Preset all GAMP defaults
if nargin <= 2 || isempty(GAMPopt)
    GAMPopt = GampOpt();
    GAMPopt.nit = 20;
    GAMPopt.tol = 1e-5;
    GAMPopt.removeMean = false;
end
GAMPopt.varNorm = false;
GAMPopt.stepWindow = 0;
GAMPopt.pvarStep = false;

%Preset all EM defaults
if nargin <= 3  || isempty(EMopt)
    EMopt.noise_var = sum(abs(Y).^2,1)/(M*101);
    EMopt.maxEMiter = 20;
    EMopt.EMtol = 1e-5;
end

%Check inputs and initializations
EMopt = check_opts(EMopt);
if EMopt.L == 1
    %Initialize BG
    [lambda, theta, phi, EMopt] = set_initsBG(EMopt, Y, A, M, N, T);
    L = 1;
else
    %Initialize GM
    [lambda, omega, theta, phi, EMopt] = set_initsGM(EMopt, Y, A, M, N, T);
    %Find number of mixture components
    L = size(theta,3);
    if (L ~= size(phi,3) || L ~= size(omega,3))
        error('There are an unequal amount of components for the active means, variances, and weights')
    end
end

histFlag = false;
if nargout >=4 
    histFlag = true;
    estHist.step = [];
    estHist.val = [];
    estHist.pass = [];
end;

muw = EMopt.noise_var;

%Initialize loop
t = 0;
stop = 0;
if ~EMopt.heavy_tailed
    tmax = ceil(3*EMopt.maxEMiter/4);
else
    tmax = EMopt.maxEMiter;
end;

%Find initialization for GAMP
if L == 1
    Xhat = ones(N,T).*theta;
else
    Xhat = ones(N,T).*sum(omega.*theta,3);
end

%If measurements are real perform real version of EMGMAMP
if isreal(Y)
    while stop == 0

        %Increment time exit loop if exceeds maximum time
        t = t + 1;
        if t >= EMopt.maxEMiter
            stop = 1;
        end
        
        %Input channel
        if L == 1
            inputEst = AwgnEstimIn(theta, phi);
        else
            inputEst = GMEstimIn(omega,theta,phi);
        end
        inputEst = SparseScaEstim(inputEst,lambda);

        %Output channel
        outputEst = AwgnEstimOut(Y, muw);   

        %Perform GAMP
        [Xhat2, Xvar, Rhat, Rvar,...
            Shat, ~, Zhat, Zvar, estHistNew] = ...
            gampEst(inputEst, outputEst, A, GAMPopt);
        
        %Output Histories if necessary
        if histFlag
            estHist.step = [estHist.step; estHistNew.step];
            estHist.val = [estHist.val; estHistNew.val];
            estHist.pass = [estHist.pass; estHistNew.pass];
        end

        if L ==1
            %Update BG parameters
            [lambda, theta, phi] = BG_update(Rhat, Rvar, lambda, theta, phi, EMopt);
        else
            %Update GM parameters
            [lambda, omega, theta, phi] = GM_update(Rhat, Rvar, lambda, omega, theta, phi, EMopt);
        end

        %Update noise variance. Include only a portion of the Zvar
        %in beginning stages of EMGMAMP because true update may make it
        %unstable.
        if EMopt.learn_noisevar
            if strcmp(EMopt.sig_dim,'joint')
                if t<tmax
                    muw = sum(sum((Y-A.mult(Xhat)).^2))/M/T;
                else
                    muw = sum(sum((Y-Zhat).^2))/M/T+sum(sum(Zvar))/M/T;
                end
            else
                if t<tmax
                    muw = sum((Y-A.mult(Xhat)).^2,1)/M;
                else
                    muw = sum((Y-Zhat).^2,1)/M+sum(Zvar,1)/M;
                end
            end
        end
        
        muw = resize(muw,M,T);

        %Calculate the change in signal estimates
        norm_change = norm(Xhat-Xhat2,'fro')^2/norm(Xhat,'fro')^2;

        %Check for estimate tolerance threshold
        if norm_change < EMopt.EMtol
            stop = 1;
        end

        %Reinitialize GAMP estimates
        Xhat = Xhat2;
        GAMPopt.xhat0 = Xhat2;
        GAMPopt.xvar0 = Xvar;
        GAMPopt.shat0 = Shat;

    end;

    %Do a final FULL EM update of noise var (psi)
    if EMopt.learn_noisevar
        if strcmp(EMopt.sig_dim,'joint')
            muw = sum(sum((Y-Zhat).^2))/M/T+sum(sum(Zvar))/M/T;
        else
            muw = sum((Y-Zhat).^2,1)/M+sum(Zvar,1)/M;
        end
    end
    muw = resize(muw,M,T);
   

    %Input channel
    if L == 1
        inputEst = AwgnEstimIn(theta, phi);
    else
        inputEst = GMEstimIn(omega,theta,phi);
    end
    inputEst = SparseScaEstim(inputEst,lambda);

    %Output channel
    outputEst = AwgnEstimOut(Y, muw);   

    %Perform GAMP
    [Xhat, Xvar, ~, ~,...
        ~, ~, ~, ~, estHistNew] = ...
        gampEst(inputEst, outputEst, A, GAMPopt);

%%
%Complex case
else
    while stop == 0

        %Increment time exit loop if exceeds maximum time
        t = t + 1;
        if t >= EMopt.maxEMiter
            stop = 1;
        end

        %Input channel
        if L == 1
            inputEst = CAwgnEstimIn(theta,phi);
        else
            inputEst = CGMEstimIn(omega,theta,phi);
        end
        inputEst = SparseScaEstim(inputEst,lambda);

        %Output channel
        outputEst = CAwgnEstimOut(Y, muw);   

        %Perform GAMP
        [Xhat2, Xvar, Rhat, Rvar,...
            Shat, ~, Zhat, Zvar, estHistNew] = ...
            gampEst(inputEst, outputEst, A, GAMPopt);
        
        if histFlag
            estHist.step = [estHist.step; estHistNew.step];
            estHist.val = [estHist.val; estHistNew.val];
            estHist.pass = [estHist.pass; estHistNew.pass];
        end

        if L ==1
            %Update BG parameters
            [lambda, theta, phi] = CBG_update(Rhat, Rvar, lambda, theta, phi, EMopt);
        else
            %Update GM parameters
            [lambda, omega, theta, phi] = CGM_update(Rhat, Rvar, lambda, omega, theta, phi, EMopt);
        end

        %Update noise variance. Include only a portion of the Zvar
        %in beginning stages of EMGMAMP because true update may make it
        %unstable.
        if EMopt.learn_noisevar
            if strcmp(EMopt.sig_dim,'joint')
                if t<tmax
                    muw = sum(sum(abs(Y-A.mult(Xhat)).^2))/M/T;
                else
                    muw = sum(sum(abs(Y-Zhat).^2))/M/T+sum(sum(Zvar))/M/T;
                end
            else
                if t<tmax
                    muw = sum(abs(Y-A.mult(Xhat)).^2,1)/M;
                else
                    muw = sum(abs(Y-Zhat).^2,1)/M+sum(Zvar,1)/M;
                end
            end
        end
        muw = resize(muw,M,T);

        %Calculate the change in signal estimates
        norm_change = norm(Xhat-Xhat2,'fro')^2/norm(Xhat,'fro')^2;

        %Check for estimate tolerance threshold
        if norm_change < EMopt.EMtol
            stop = 1;
        end

        %Reinitialize GAMP estimates
        Xhat = Xhat2;
        GAMPopt.xhat0 = Xhat2;
        GAMPopt.xvar0 = Xvar;
        GAMPopt.shat0 = Shat;

    end;

        %Do a final FULL EM update of noise var (psi)
        if EMopt.learn_noisevar
            if strcmp(EMopt.sig_dim,'joint')
                muw = sum(sum(abs(Y-Zhat).^2))/M/T+sum(sum(Zvar))/M/T;
            else
                muw = sum(abs(Y-Zhat).^2,1)/M+sum(Zvar,1)/M;
            end
        end
    muw = resize(muw,M,T);

    %Input channel
    if L == 1
        inputEst = CAwgnEstimIn(theta,phi);
    else
        inputEst = CGMEstimIn(omega,theta,phi);
    end
    inputEst = SparseScaEstim(inputEst,lambda);

    %Output channel
    outputEst = CAwgnEstimOut(Y, muw);   

    %Perform GAMP
    [Xhat, Xvar, ~, ~,...
        ~, ~, ~, ~, estHistNew] = ...
        gampEst(inputEst, outputEst, A, GAMPopt);
    
end

%Output Final Histories
if histFlag
    estHist.step = [estHist.step; estHistNew.step];
    estHist.val = [estHist.val; estHistNew.val];
    estHist.pass = [estHist.pass; estHistNew.pass];
end

%Output final parameters
%BG parameters
if L ==1
    if strcmp(EMopt.sig_dim,'joint')
        param.lambda = lambda(1,1);
        param.active_mean = theta(1,1);
        param.active_var = phi(1,1);
        param.noise_var = muw(1,1);
    else
        param.lambda = lambda(1,:);
        param.active_mean = theta(1,:);
        param.active_var = phi(1,:);
        param.noise_var = muw(1,:);
    end
%GM parameters
else
    if strcmp(EMopt.sig_dim,'joint')
        param.lambda = lambda(1,1);
        param.active_weights = reshape(omega(1,1,:),1,L)';
        param.active_mean = reshape(theta(1,1,:),1,L).';
        param.active_var = reshape(phi(1,1,:),1,L)';
        param.noise_var = muw(1,1);
    else
        param.lambda = lambda(1,:);
        param.active_weights = reshape(omega(1,:,:),T,L)';
        param.active_mean = reshape(theta(1,:,:),T,L).';
        param.active_var = reshape(phi(1,:,:),T,L)';
        param.noise_var = muw(1,:);
    end
end    

return;

function [xhatFinal, xvarFinal, rhatFinal, rvarFinal,...
    shatFinal, svarFinal, zhatFinal,zvarFinal,estHist] = ...
    gampEst(scaEstIn, scaEstOut, A, opt)
% gampEstRobust:  G-AMP estimation algorithm (robust implementation)
%
% The G-AMP estimation algorithm is intended for the estimation of a
% random vector x observed through an observation y from the Markov chain
%
%   x -> z = A*x -> y,
%
% where the components of x are independent and the mapping z -> y is
% separable.  The function takes three arguments:
%
% scaEstIn:  An input estimator derived from the EstimIn class
%    based on the input distribution p_X(x_j).
% scaEstOut:  An output estimator derived from the EstimOut class
%    based on the output distribution p_{Y|Z}(y_i|z_i).
% A:  Either a matrix or a linear operator defined by the LinTrans class.
% opt:  A set of options of the class GampOpt.
%
% xhatFinal: The final estimate of the vector x.
% rhatFinal, rvarFinal:  Final estimate for rhat and rvar (see entries in
%     estHist below).
%
% estHist:  History of the estimator per iteration.  Note that this history
%    includes steps that may have been aborted.
% estHist.rhat, estHist.rvar:  When the estimator functions implement
%    the sum-product algorithm, the conditional distribution of x(j)
%    given the vector y is  approximately
%       p(x(j)|y) = C*p(x(j))*exp( -(rhat(j)-x(j))^2/(2*rvar(j) )
%    where C is a normalization constant.
% estHist.val:  The history of the value function.
% estHist.xhat:  The history of the estimates.
% estHist.xvar:  The history of the estimates.
% estHist.scaleFac:  The history of the variance-normalization scaling factor.

% Get options
if (nargin < 4)
    opt = GampOpt();
elseif (isempty(opt))
    opt = GampOpt();
end
nit     = opt.nit;              % number of iterations
step    = opt.step;             % step size
stepMin = opt.stepMin;          % minimum step size
stepMax = opt.stepMax;          % maximum step size
stepIncr = opt.stepIncr;	% step increase
stepDecr = opt.stepDecr;	% step decrease
adaptStep = opt.adaptStep;      % adaptive step size
stepWindow = opt.stepWindow;    % step size check window size
bbStep = opt.bbStep;            % Barzilai Borwein step size
verbose = opt.verbose;          % Print results in each iteration
tol = opt.tol;                  % Convergence tolerance
stepTol = opt.stepTol;          % minimum allowed step size
compVal = true;
saveHist = (nargout >= 9);
pvarStep = opt.pvarStep;        % incldue step size in pvar/pvarBar
varNorm = opt.varNorm;          % normalize variances
scaleFac = opt.scaleFac;	% initial variance normalization 
%Make variances uniform by averaging. This options does not currently work
%with mean removal!
uniformVariance = opt.uniformVariance;

% If A is a double matrix, replace by an operator
if isa(A, 'double')
    A = MatrixLinTrans(A);
end

% Get dimensions
[m,n] = A.size();
S = scaEstOut.numColumns();


% Get default initialization values 
[xhat0,xvar0] = scaEstIn.estimInit();
valIn = -inf;       % Initialize step-sizing quantity

%Replace default initialization with warm start values if provided 
if ~isempty(opt.xhat0)
    xhat0 = opt.xhat0;
end
if ~isempty(opt.xvar0)
    xvar0 = opt.xvar0;
end

% For a scalar output, the same distribution is applied to all components
if (length(xhat0) == 1)
    xhat0 = repmat(xhat0,n,S);
end
if (length(xvar0) == 1)
    xvar0 = repmat(xvar0,n,S);
end

% If the mean removal option is set, create an augmented linear
% transform class with the mean removed.  (See the LinTransDemean
% class for more details.)  Then initialize xhat and xvar
if (opt.removeMean)
    Ad = LinTransDemean(A,xhat0);
    xvar0Mean = mean(xvar0);
    xhat = zeros(n+1,1);
    xvar = [xvar0; xvar0Mean];
    rhat = xhat;
    rvar = 100*xvar;
else
    Ad = A;
    xhat = xhat0;
    xvar = xvar0;
    rhat = nan(size(xhat));
    rvar = nan(size(xhat));
end
[m1,n1] = Ad.size();

% Continue with initialization 
shat = zeros(m1,S);	% default value is zero
svar = nan(m1,S);	% will test for NaN later
xhatPrev = nan(n1,S);	% will test for NaN later

%Replace default initialization with warm start values if provided 
if ~isempty(opt.shat0)
    shat = opt.shat0*scaleFac;	% variance normalization included
end
if ~isempty(opt.svar0)
    svar = opt.svar0*scaleFac; 	% variance normalization included
end
if ~isempty(opt.xhatPrev0)
    xhatPrev = opt.xhatPrev0;
end

% Declare variables
%zhat0 = Ad.mult(xhat);
zhat0 = nan(m1,S);
zvar0 = nan(size(zhat0));
xhatFinal  = nan(n1,S);
val = nan(nit,1);
valOpt = [];
if (saveHist)
    estHist.xhat = nan(n*S,nit);
    estHist.rhat = nan(n*S,nit);
    estHist.phat = nan(m*S,nit);
    if ~uniformVariance
        estHist.xvar = nan(n*S,nit);
        estHist.rvar = nan(n*S,nit);
        estHist.pvar = nan(m*S,nit);
    else
        estHist.xvar = nan(1,nit);
        estHist.rvar = nan(1,nit);
        estHist.pvar = nan(1,nit);
    end
    estHist.val = nan(nit,1);
    estHist.step = nan(nit,1);
    estHist.pass = nan(nit,1);
    estHist.scaleFac = nan(nit,1);
end


% Check for the presence of two methods within the LinTrans and EstimIn
% objects and set flags accordingly
MtxUncertaintyFlag = ismethod(Ad,'includeMatrixUncertainty');
MsgUpdateFlag = ismethod(scaEstIn, 'msgUpdate');

%For uniform variances, we need a scaled version of the squared Frobenius
%norm of A
if uniformVariance
    %Compute squared Frobenius norm and normalize by number of elements
    mA2 = sum(Ad.multSq(ones(n,1)))/m/n;

    %If Ad.multSq returned NaN, then use a different technique
    if isnan(mA2),
      numProbe = 10;
      mA2 = 0;
      for i=1:numProbe
        mA2 = mA2 + (norm(A.mult(randn(n,1)),'fro')^2)/m/n/numProbe;
      end;
    end 
end

% Convert relative pvarMin and xvarMin to actual minimum quantities
if uniformVariance
    xvar = mean(xvar(:));
    rvar = mean(rvar(:));
    svar = mean(svar(:));
    pvarMin = opt.pvarMin*n*mA2*mean(xvar0(:)); %Asq(1:m,1:n)*xvar0; ---
else 
    pvarMin = opt.pvarMin*(A.multSq(xvar0)); %Asq(1:m,1:n)*xvar0; ---
end
pvarMin(pvarMin == 0) = opt.pvarMin; %Protect from zero valued pvarMin
xvarMin = opt.xvarMin*xvar;

%If computing BB steps, compute column norms for use in scaling
if bbStep
    columnNorms = A.multSqTr(ones(m,1)).^0.5;
    columnNorms = repmat(columnNorms,1,S);
end

%Control variable to end the iterations
stop = false;
it = 0;

%Init to dummy value
step1 = 1;		% over-rides step in first half of first iteration
pvarOpt = 0;		% value is inconsequential when step1=1
pvarBarOpt = 0;		% value is inconsequential when step1=1


% Main iteration loop
while ~stop
    
    % Iteration count
    it = it + 1;
    
    % Check for final iteration
    if it >= nit
        stop = true;
    end
    
    % Output linear step with no A uncertainty
    if ~uniformVariance
        pvarBar = Ad.multSq(xvar);
    else
        %Verify that xvar is scalar
        mxvar = norm(xvar(:),1)/numel(xvar);
        %Scale with normalized squared Frobenius norm of A, times number of
        %unknowns n
        pvarBar = n*mA2*mxvar;
    end
    
    %Incorporate A uncertainty
    if MtxUncertaintyFlag
        pvar = Ad.includeMatrixUncertainty(pvarBar,xhat,xvar);
    else
        pvar = pvarBar;
    end
    
    %Step in pvar
    if pvarStep
        pvar = step1*pvar + (1-step1)*pvarOpt;
        pvarBar = step1*pvarBar + (1-step1)*pvarBarOpt;
    end
    
    %Update zhat
    zhat = Ad.mult(xhat);
    
    % Compute log likelihood at the output and add it the total negative
    % K-L distance at the input.
    if (compVal)
        if ~uniformVariance
            valOut = sum(sum(scaEstOut.logLike(zhat(1:m,:),pvar(1:m,:))));
        else
            valOut = sum(sum(scaEstOut.logLike(zhat(1:m,:),pvar)));
        end
        val(it) = valOut + valIn;
    end
    
    % Determine if candidate passed
    if (it > 1)
        
        %Check against worst value in last stepWindow good steps
        stopInd = length(valOpt);
        startInd = max(1,stopInd - stepWindow);
        
        %Check the step
        pass = (val(it) > min(valOpt(startInd:stopInd))) ||...
            ~adaptStep || (step <= stepMin);
        
    else
        pass = true;
    end
    
    %Save the step size and pass result if history requested
    if saveHist
        estHist.step(it) = step;
        estHist.pass(it) = pass;
    end
    
    % If pass, set the optimal values and compute a new target shat and
    % snew.
    if (pass)
        
        % Set new optimal values
        shatOpt = shat;
        svarOpt = svar;
        xhatPrevOpt = xhatPrev;
        xhatOpt = xhat;
        pvarOpt = pvar;
        pvarBarOpt = pvarBar;
        
        %Set zhat outputs
        zhatFinal = zhat0;
        zvarFinal = zvar0;
        
        %We keep a record of only the succesful step valOpt values
        valOpt = [valOpt val(it)]; %#ok<AGROW>
        
        %Store previous optimal solution
        xhatPrevFinal = xhatFinal;
        
        % Save current optimal solution
        if (opt.removeMean)
            [xhatFinal, rhatFinal, rvarFinal] = Ad.getEst(xhat,rhat,rvar*scaleFac);
            xvarFinal = xvar(1:n);
        else
            xhatFinal = xhatOpt;
            xvarFinal = xvar;
            rhatFinal = rhat;
            rvarFinal = rvar*scaleFac;	% report unscaled version
        end
        
        %Set shat and svar outputs
        shatFinal = shatOpt/scaleFac;	% report unscaled version
        svarFinal = svarOpt/scaleFac;	% report unscaled version
        
        % Check for convergence
        if (it > 1) && ...
                (norm(xhatPrevFinal - xhatFinal) / norm(xhatFinal) < tol)
            stop = true;
        end
        
        % Continued output step
        phat = zhat - (pvarBar/scaleFac).*shat; %uses pvarBar rather than pvar
        if ~uniformVariance
            pvar(1:m,:) = max(pvar(1:m,:), pvarMin);
        else
            pvar = max(pvar,pvarMin);
        end
        
        %Set scaleFac to mean of pvar if normalization is on. 
	%Else scaleFac remains at the initialized value of 1 and has no effect
        if varNorm
            scaleFac = mean(pvar(:));
        end
        
        % Output nonlinear step -- regular components
        if ~uniformVariance
            shatNew = zeros(m1,S);
            svarNew = zeros(m1,S);
            I = (1:m)';
            [zhat0,zvar0] = scaEstOut.estim(phat(I,:),pvar(I,:));
            shatNew(I,:) = (scaleFac./pvar(I,:)).*(zhat0-phat(I,:));
            svarNew(I,:) = (scaleFac./pvar(I,:)).*(1-zvar0./pvar(I,:));
        else
            [zhat0,zvar0] = scaEstOut.estim(phat,pvar);
            shatNew = (scaleFac./pvar).*(zhat0-phat);
            svarNew = (scaleFac./pvar).*(1-zvar0./pvar);
        end
        
        % Output nonlinear step for the demean output component
        if (opt.removeMean)
            shatNew(m+1) = (-phat(m+1))/(pvar(m+1));
            svarNew(m+1) = 1/pvar(m+1);
        end
        
        %Scalar Variance
        if uniformVariance
            svarNew = mean(svarNew(:));
        end
        
        
        %Compute new BB Step size if requested
        if bbStep && it > 2
            %Compute previous step direction/size weighted with column
            %norms
            sBB = (xhatOpt(1:n,:) - xhatPrevOpt(1:n,:));
            
            %Compute new step size using columnNorms weighting
            %Select the smallest step over all the columns for a matrix
            %valud signal
            values = sum(abs(sBB .* columnNorms).^2,1) ./...
                sum(abs(A.mult(sBB).^2),1);
            step = min(values);
        end
        
        %Enforce step size bounds
        step = min([stepIncr*max([step stepMin]) stepMax]);
        
    else % not pass
        % Decrease step size
        step = max(stepMin, stepDecr*step);
        
        %Check for minimum step size
        if step < stepTol
            stop = true;
        end
    end % if pass
    
    % Save results
    if (saveHist)
        estHist.xhat(:,it) = xhatFinal(:);
        estHist.rhat(:,it) = rhatFinal(:);
        estHist.phat(:,it) = reshape(phat(1:m,:),[],1);
        if ~uniformVariance
            estHist.pvar(:,it) = reshape(pvar(1:m,:),[],1);
            estHist.xvar(:,it) = xvarFinal(:);
            estHist.rvar(:,it) = rvarFinal(:);
        else
            estHist.pvar(it) = pvar;
            estHist.xvar(it) = xvarFinal;
            estHist.rvar(it) = rvarFinal;
        end
        estHist.val(it) = val(it);
        estHist.scaleFac(it) = scaleFac;
    end
    
    % Create new candidate shat
    step1 = step;
    if (it==1)
      if any(isnan(svarOpt)), % if user didn't specify svar0 
        svarOpt = svarNew; 	% equivalent to making step1=1
      end;
      if any(isnan(xhatPrevOpt)), % if user didn't specify xhatPrev0
        xhatPrevOpt = xhatOpt; % equivalent to making step1=1
      end;
    end;
    shat = (1-step1)*shatOpt + step1*shatNew;
    svar = (1-step1)*svarOpt + step1*svarNew;
    xhatPrev = (1-step1)*xhatPrevOpt + step1*xhatOpt;
    
    % Print results
    if (verbose)
        fprintf(1,'it=%3d value=%12.4e step=%f\n', it, val(it), step);
    end
    
    % Input linear step
    if ~uniformVariance
        rvar = 1./Ad.multSqTr(svar);   % rvar = 1./((A.^2)*svar)
    else
        rvar = 1/(mA2*svar*m); %use Frobenius norm with scalar variance
    end
    
    %Rest of input linear step
    rhat = xhatPrev + rvar.*(Ad.multTr(shat)); % rhat = xhatPrev + rvar.*(A*shat)
    rvar = max(rvar, xvarMin);
    
    % Send messages to input estimation function.
    if MsgUpdateFlag
        valMsg = scaEstIn.msgUpdate(it, rhat, rvar);
    else
        valMsg = 0;
    end
    
    % Input nonlinear step
    if (opt.removeMean)
        
        % Regular components
        I = (1:n)';
        [xhat(I),xvar(I),valIn] = scaEstIn.estim(rhat(I)+xhat0, rvar(I));
        xhat(I) = xhat(I) - xhat0;
        
        % Mean component
        xhat(n+1) = xvar0Mean/(xvar0Mean+rvar(n+1))*rhat(n+1);
        xvar(n+1) = xvar0Mean*rvar(n+1)/(xvar0Mean+rvar(n+1));
        
    else
        % Call input scalar estimator
        [xhat,xvar,valIn] = scaEstIn.estim(rhat, rvar*scaleFac);
        
    end
    
    %Scalar variances
    if uniformVariance
        xvar = mean(xvar(:));
    end
    
    
    valIn = sum( valIn(:) ) + valMsg;
    
end

%Trim the outputs if early termination occurred
if saveHist && (it < nit)
    estHist.xhat = estHist.xhat(:,1:it);
    estHist.rhat = estHist.rhat(:,1:it);
    estHist.phat = estHist.phat(:,1:it);
    estHist.val = estHist.val(1:it);
    estHist.scaleFac = estHist.scaleFac(1:it);
    estHist.step = estHist.step(1:it);
    estHist.pass = estHist.pass(1:it);
    if ~uniformVariance
        estHist.rvar = estHist.rvar(:,1:it);
        estHist.pvar = estHist.pvar(:,1:it);
        estHist.xvar = estHist.xvar(:,1:it);
    else
        estHist.rvar = estHist.rvar(1:it);
        estHist.pvar = estHist.pvar(1:it);
        estHist.xvar = estHist.xvar(1:it);
    end
end


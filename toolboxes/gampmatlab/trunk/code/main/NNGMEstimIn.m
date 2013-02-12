classdef NNGMEstimIn < EstimIn
    % NNGMEstimIn:  Non-negativeGaussian Mixture 
    % scalar input estimation function
    
    properties 
        % Prior mean and variance
        omega; % Weights
        theta;  % Means 
        phi;   % Variances 
        const; %normalizing constant
    end
    
    methods
        % Constructor
        function obj = NNGMEstimIn(omega, theta, phi)
            obj = obj@EstimIn;
            obj.omega = omega;
            obj.theta = theta;
            obj.phi = phi;
            obj.const = max(erfc(-theta./sqrt(phi)/sqrt(2)),1e-300)/2;
        end

        % Prior mean and variance
        function [mean0, var0, valInit] = estimInit(obj)
            %Define necessary constants
            alpha = -obj.theta./sqrt(obj.phi);
            inv_mill = sqrt(2/pi)./erfcx(alpha/sqrt(2));
            %Find mean and variance of the prior
            mean0 = sum(obj.omega.*erfc(alpha/sqrt(2))/2.*...
                (obj.theta + sqrt(obj.phi).*inv_mill),3)./...
                sum(obj.omega.*erfc(alpha/sqrt(2))/2,3);
            var0  = sum(obj.omega.*erfc(alpha/sqrt(2))/2.*...
                (obj.phi.*(1-inv_mill.*(inv_mill - alpha))...
                + (obj.theta + sqrt(obj.phi).*inv_mill).^2),3)./...
                sum(obj.omega.*erfc(alpha/sqrt(2))/2,3) - mean0.^2;
            valInit = 0;
        end


        function [Uhat, Uvar, NKL] = estim(obj, Rhat, Rvar)
            
            Rhat = real(Rhat);
            
            %Get the number of mixture components
            L = size(obj.theta,3);
            %Find the signal dimension
            [M, N] = size(Rhat);
            
            if L == 1
                
                beta = obj.phi + Rvar + eps;
                gamma = (Rvar.*obj.theta+obj.phi.*Rhat)./beta;
                nu = obj.phi.*Rvar./beta;
                eta = -gamma./sqrt(nu);

                cdf_comp = max(erfc(eta/sqrt(2)),1e-300)/2;
                inv_mill = sqrt(2/pi)./erfcx(eta/sqrt(2));

                % Compute MMSE estimates
                Uhat = gamma + sqrt(nu).*inv_mill;
                Uvar = nu.*(1-inv_mill.*(inv_mill - eta))...
                    + (gamma + sqrt(nu).*inv_mill).^2 - Uhat.^2;
                
                % Compute the negative KL divergence            
                if (nargout >= 3)                            
                    NKL0 = 0.5*log(Rvar./beta./(obj.const).^2)...
		    	-(Rhat-obj.theta).^2./2./beta...
			+(Uvar +(Uhat - Rhat).^2)./(2*Rvar);
		    NKL = NKL0 + log(cdf_comp);
                    %Find indices that could cause numerical issues in 
                    %erfc computation and use erfc(x) \approx 
		    %(0.3480242*t-0.0958798*t^2+0.7478556*t^3)*exp(-x^2)
		    %for t=1/(1+.47047*x), from Abramowitz and Stegun
                    I = find(eta > 10);
		    tI = 1./(1+0.47047*eta(I)/sqrt(2));
		    NKL(I) = NKL0(I) -0.5*eta(I).^2 + log(0.5)...
		    	+log(0.3480242*tI-0.0958798*tI.^2+0.7478556*tI.^3);
                end
            
            else
            
                obj.omega = obj.omega ./ repmat(sum(obj.omega, 3), [1, 1, L]);


                %Preallocate storage
                eta = zeros(M,N,L); gamma = zeros(M,N,L); 
                nu = zeros(M,N,L); alpha = zeros(M,N,L);

                %Run through mixture components
                for i = 1:L
                    eta(:,:,i) = obj.phi(:,:,i) + Rvar + eps;
                    gamma(:,:,i) = (Rvar.*obj.theta(:,:,i)+ ...
                        obj.phi(:,:,i).*Rhat)./eta(:,:,i);
                    alpha(:,:,i) = (Rhat-obj.theta(:,:,i)).^2./eta(:,:,i);
                    nu(:,:,i) = obj.phi(:,:,i).*Rvar./eta(:,:,i);
                end;   

                %compute terms inside of cdf/pdf components
                beta = -gamma./sqrt(nu);
                cdf_comp = max(erfc(beta/sqrt(2)),1e-300);
                inv_mill = sqrt(2/pi)./erfcx(beta/sqrt(2));

                lik = zeros(M,N,L);
                for i = 1:L
                    lik = lik + cdf_comp./repmat(cdf_comp(:,:,i),[1 1 L]).*...
                        repmat(obj.omega(:,:,i),[1 1 L])./obj.omega...
                        .*sqrt(eta./repmat(eta(:,:,i),[1 1 L]))...
                        .*exp((alpha - repmat(alpha(:,:,i),[1 1 L]))/2);
                end

                %Compute MMSE quantities
                Uhat = sum((gamma + sqrt(nu).*inv_mill)./lik,3);
                Uvar = sum((nu.*(1-inv_mill.*(inv_mill - beta))...
                    + (gamma + sqrt(nu).*inv_mill).^2)./lik,3) - Uhat.^2;


                % Compute the negative KL divergence            
                if (nargout >= 3)                            
                    NKL = log(sum(obj.omega.*exp(-alpha/2)./sqrt(beta).*cdf_comp/2,3)...
                        ./sum(obj.omega.*obj.const,3)...
                        .*sqrt(Rvar))+(Uvar +(Uhat - Rhat).^2)./(2*Rvar);
                end
            
            end

        end
        
        % Generate random samples
        function x = genRand(obj, nx, nt)
                L = size(obj.theta,3);
                dummy2 = cumsum(obj.omega,3);
                dummy = zeros(nx, nt, L+1);
                dummy(:,:,2:end) = dummy2;
                dummy2 = rand(nx,nt);
                
                x = zeros(nx,nt);
                
                
                for i = 1:L
                    temp = erfc(-obj.theta(:,:,i)./sqrt(2.*obj.phi(:,:,i)))/2;
                    temp2 = (dummy2>=dummy(:,:,i) & dummy2<dummy(:,:,i+1));
                    temp = sqrt(2.*obj.phi(:,:,i)).*erfinv(2*(1- temp...
                        + rand(nx,nt,1).*temp)-1)+obj.theta(:,:,i);
                    x(temp2) = temp(temp2);
                end

        end
        
        % Computes the likelihood p(y) for y = x + v, v = N(0,Yvar)
        function py = plikey(obj,Y,Yvar)
            
            Y = real(Y);
            
            L = size(obj.omega,3);
            [M, T] = size(Y);
            lik = zeros(M,T,L);
            gamma = zeros(M,T,L); nu = zeros(M,T,L);
            
            for i = 1:L
                dummy = obj.phi(:,:,i) + Yvar;
                gamma(:,:,i) = (Yvar.*obj.theta(:,:,i)+obj.phi(:,:,i).*Y)./dummy;
                nu(:,:,i) = obj.phi(:,:,i).*Yvar./dummy;
                lik(:,:,i) = obj.omega(:,:,i).*exp(-1./(2*(obj.phi(:,:,i)+Yvar+eps))...
                    .*(Y-obj.theta(:,:,i)).^2)./sqrt(2*pi*(obj.phi(:,:,i)+Yvar+eps));
            end
            
            lik = sum(lik.*erfc(-gamma./sqrt(nu)/sqrt(2))/2,3)...
                ./sum(obj.omega.*obj.const,3);
            
            lik(isnan(lik)) = 0.999;
            py = sum(lik,3);
        end

    end
    
end


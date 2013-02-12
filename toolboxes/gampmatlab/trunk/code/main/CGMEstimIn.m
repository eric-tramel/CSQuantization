classdef CGMEstimIn < EstimIn
    % CGMEstimIn:  Complex Gaussian Mixture scalar input estimation function
    
    properties 
        % Prior mean and variance
        omega; % Weights
        theta;  % Means 
        phi;   % Variances 
    end
    
    methods
        % Constructor
        function obj = CGMEstimIn(omega, theta, phi)
            obj = obj@EstimIn;
            obj.omega = omega;
            obj.theta = theta;
            obj.phi = phi;            
        end

        % Prior mean and variance
        function [mean0, var0, valInit] = estimInit(obj)
            mean0 = sum(obj.theta.*obj.omega,3);
            var0  = sum(obj.omega .* (obj.phi + ...
                abs(obj.theta).^2), 3) - abs(mean0).^2;
            valInit = 0;
        end 


        function [Uhat, Uvar, NKL] = estim(obj, Rhat, Rvar)

            %Get the number of mixture components
            L = size(obj.omega,3);
            obj.omega = obj.omega ./ repmat(sum(obj.omega, 3), [1, 1, L]);
            
            %Grab the signal dimension
            [N, T] = size(Rhat);
            
            %Preallocate storage
            gamma = zeros(N,T,L); alpha = zeros(N,T,L);
            beta = zeros(N,T,L); nu = zeros(N,T,L);

            for i = 1:L
               beta(:,:,i) = obj.phi(:,:,i) + Rvar + eps;
               alpha(:,:,i) = abs(Rhat-obj.theta(:,:,i)).^2./beta(:,:,i);
               gamma(:,:,i) = (Rhat.*obj.phi(:,:,i) + obj.theta(:,:,i).*Rvar)./beta(:,:,i);
               nu(:,:,i) = Rvar.*obj.phi(:,:,i)./beta(:,:,i);
            end

            lik = zeros(N,T,L);
            for i = 1:L
                lik = lik + repmat(obj.omega(:,:,i),[1 1 L])./obj.omega...
                    .*beta./repmat(beta(:,:,i),[1 1 L])...
                    .*exp((alpha-repmat(alpha(:,:,i),[1 1 L])));
            end

            Uhat = sum(gamma./lik,3);
            Uvar = sum((nu + abs(gamma).^2)./lik,3) - abs(Uhat).^2;
            

            % Compute the negative KL divergence            
            if (nargout >= 3)                            
                NKL = log(sum(obj.omega.*exp(-alpha)./beta,3)...
                    .*Rvar)+(Uvar + abs(Uhat - Rhat).^2)./(2*Rvar);
            end

        end
        
        % Generate random samples
        function x = genRand(obj, nx)            
                L = length(obj.theta);
                dummy = [0 cumsum(obj.omega)];
                dummy2 = rand(nx,1);

                dummy3 = zeros(nx,1,L);
                for i = 1:L
                    dummy3(:,i) = (dummy2>=dummy(i) & dummy2<dummy(i+1))...
                        .*(obj.theta(i) + sqrt(obj.phi(i)/2)...
                        .*(randn(nx,1)+1i*randn(nx,1)));
                end

                x = sum(dummy3,2);
        end
        
        % Computes the likelihood p(y) for y = x + v, v = N(0,Yvar)
        function py = plikey(obj,Y,Yvar)
            
            L = size(obj.omega,3);
            [M, T] = size(Y);
            lik = zeros(M,T,L);
            
            for i = 1:L
                lik(:,:,i) = obj.omega(:,:,i).*exp(-1./((obj.phi(:,:,i)+Yvar))...
                    .*abs(Y-obj.theta(:,:,i)).^2)./(2*pi*(obj.phi(:,:,i)+Yvar));
            end
            
            py = sum(lik,3);
        end
        
        % Computes the log-likelihood, log p(Y(i,j)), for Y = X + V, where 
        % p(X(i,j)) = sum_k omega(i,j,k)*CN(theta(i,j,k), phi(i,j,k)) and 
        % p(V(i,j)) = CN(0, Yvar(i,j))
        function logpy = loglikey(obj, Y, Yvar)
            logpy = log(obj.plikey(Y, Yvar));
        end

    end
    
end


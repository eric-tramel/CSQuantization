classdef MriLinTrans < LinTrans
    % Linear transformation for a MRI
    %
    % This implements the function z = A(x) where A is the block matrix
    %
    %   A = [F; W; TV]
    %
    % where F = Fourier transform * coil sensitivities
    %       W = Wavelet
    %       TV = total variation
    %
    % Right now the code only implements the Fourier part
    
    properties %(Access = private)
        % Parameter structure used by the Stanford code
        param;
        
        % Dimensions
        nrow, ncol;     % number of rows and columns
        
        % Fourier parameters
        ncoil;      % number of coils
        nsamp;      % number of rows after sub-sampling
                    % note that each col is fully sampled
        sampMask;   % sampMask(j)=1 if row j is sampled
        Isamp;      % indices of rows that are sampled
        noutF;      % number of outputs of the Fourier
               
        noutwav;    % number of outputs of the Wavelet
        nouttotalv; % number of outputs of the Total Variation
       
             
    end
    
    methods 
        % Constructor based on the parameter structured created
        % by the Stanford code.
        function obj = MriLinTrans(param)
            
            %%%%%%%%%%%%%%%%%%% Fourier part
            % Base class
            obj = obj@LinTrans;
            
            % Save parameter
            obj.param = param;
            
            % Get dimensions
            obj.nrow = size(param.y,1); 
            obj.ncol = size(param.y,2);
            obj.ncoil = size(param.y,3);
                    
            paramElement = getElement(param.E);
            obj.sampMask = paramElement{2}(:,1);
            
            obj.Isamp = find(obj.sampMask);
            obj.nsamp = length(obj.Isamp);            
            obj.noutF = obj.nsamp * obj.ncoil * obj.ncol;
            %%%%%%%%%%%%%%%%%%%
            
            % Wavelet 
            obj.noutwav = obj.nrow * obj.ncol;              
            
            % TV.Note the multiplication by 2 for the row and col
            % difference
            obj.nouttotalv = 2 * obj.nrow * obj.ncol;    
                     
        end       

        % Size
        function [nout,nin] = size(obj)
           nin = obj.nrow*obj.ncol;
           nout = obj.noutF + obj.noutwav +  obj.nouttotalv;
        end
        
        % extract the components sizes
        function [noutF, noutwav, nouttotalv] = getCompSize(obj)
            noutF = obj.noutF;
            noutwav = obj.noutwav;
            nouttotalv = obj.nouttotalv;
        end
        

        % Matrix multiply:  z = A*x
        function [z] = mult(obj,x)
            
            % Reshape x to be square
            xsq = reshape(x, obj.nrow, obj.ncol);
            
            % Multiply using parmeter 
            zsqF = obj.param.E * xsq;
            
            % Subsample
            zsqF = zsqF(obj.Isamp,:,:);
            
            % Return to vector
            zF = zsqF(:);
            
            % Wavelet
            zsqW = obj.param.W*xsq;            
            zW = zsqW(:);
            
            % TV
            zsqTV = obj.param.TV*xsq;
            zTV = zsqTV(:);
            
            % Place into single column
            z = [zF; zW; zTV];
                  
          
        end
        
        % Multiply by square:  pvar = abs(Ad).^2*xvar
        function [pvar] = multSq(obj, xvar)
            
            xvarsq = reshape(xvar,obj.nrow,obj.ncol);
            
            
            %%%%%%%%%%%%%%%%%%%              
            FourierElem = getElement(obj.param.E);
            FourierSenElem = FourierElem{3};

          
            FourierSenPO = abs(FourierSenElem).^2;
            pvarsqZ = zeros(obj.nrow,obj.ncol,obj.ncoil);
            for ch=1:size(FourierSenElem,3)
                pvarsqZ(:,:,ch)=(xvarsq.*FourierSenPO(:,:,ch));
            end 
            pvarsqZ = pvarsqZ(obj.Isamp,:,:);
            pvarZ = pvarsqZ(:);
            %%%%%%%%%%%%%%%%%%%
            
            
            %%%%%%%%%%%%%%%%%%%
            pvarW = sum(xvarsq(:))/(obj.nrow*obj.ncol)*ones(obj.noutwav,1);
            %%%%%%%%%%%%%%%%%%%
            
            
            %%%%%%%%%%%%%%%%%%%
            xvarrsh1 = xvarsq([2:end],:) + xvarsq([1:end-1],:);
            pvarsqTV1 = cat(1,xvarrsh1,zeros(1,obj.ncol));
            
            xvarrsh2 = xvarsq(:,[2:end]) + xvarsq(:,[1:end-1]);
            pvarsqTV2 = cat(2,xvarrsh2,zeros(obj.nrow,1));
            
            pvarsqTV = cat(3, pvarsqTV1, pvarsqTV2);
            pvarTV = pvarsqTV(:);
            %%%%%%%%%%%%%%%%%%%
            
            pvar = [pvarZ;pvarW;pvarTV];
            
        end
        

        % Matrix multiply transpose:  x = A'*s
        function [x] = multTr(obj,s)
            
            sF = s(1:obj.noutF);
            sW = s(obj.noutF+1:  obj.noutF + obj.noutwav );
            sTV = s(obj.noutF + obj.noutwav + 1 : obj.noutF + obj.noutwav + obj.nouttotalv);                     
%             
%             sF = s(1:obj.noutF);
            %%%%%%%%%%%%%%%%%%%
            ssqZ = zeros(obj.nrow,obj.ncol,obj.ncoil);
            srshZ = reshape(sF, obj.nsamp, obj.ncol, obj.ncoil);
            ssqZ(obj.Isamp,:,:) = srshZ;
            
            % Perform transform
            xsqZ = obj.param.E' * ssqZ;
            
            % Reshape
            xZ = xsqZ(:);
            %%%%%%%%%%%%%%%%%%%
            
            
            %%%%%%%%%%%%%%%%%%%
            ssqW = reshape(sW,obj.nrow,obj.ncol);
            
                        
            % Perform transform
            xsqW = obj.param.W' * ssqW;
            
            % Reshape
            xW = xsqW(:);
            %%%%%%%%%%%%%%%%%%%
            
            
            %%%%%%%%%%%%%%%%%%%
            ssqTV = reshape(sTV,obj.nrow,obj.ncol,2);

            % Perform transform
            xsqTV = obj.param.TV' * ssqTV;
            
            % Reshape
            xTV = xsqTV(:);
            %%%%%%%%%%%%%%%%%%%
            
            x = [xZ;xW;xTV];
%             x = xZ; 
        end
                    
 
        % Matrix multiply with componentwise square transpose:  
        %   rvar = (Ad.^2)'*svar
        function [rvar] = multSqTr(obj,svar)
            
            svarF = svar(1:obj.noutF);
            svarW = svar(obj.noutF+1:  obj.noutF + obj.noutwav );
            svarTV = svar(obj.noutF + obj.noutwav + 1 : obj.noutF + obj.noutwav + obj.nouttotalv);
%                         
%             svarF = svar(1:obj.noutF);
           
            %%%%%%%%%%%%%%%%%%%
            svarsqZ = zeros(obj.nrow,obj.ncol,obj.ncoil);
            svarrshZ = reshape(svarF, obj.nsamp, obj.ncol, obj.ncoil);
            svarsqZ(obj.Isamp,:,:) = svarrshZ;
            
            FourierElem = getElement(obj.param.E);
            FourierSenElem = FourierElem{3};
                     
            FourierSenPO = abs(FourierSenElem).^2;
            rvarsqZ = zeros(obj.nrow,obj.ncol);
            for ch=1:size(FourierSenElem,3)    
                rvarsqZ = rvarsqZ + FourierSenPO(:,:,ch)*sum(sum(svarsqZ(:,:,ch)))/(obj.ncol * obj.nrow);
            end 
            rvarZ = rvarsqZ(:);
            %%%%%%%%%%%%%%%%%%%
            
            
               
            %%%%%%%%%%%%%%%%%%%
            svarsqW = reshape(svarW,obj.nrow,obj.ncol);
            
            %%%%%%%%%%%%%%%%%%%             
            rvarW = sum(svarsqW(:))/(obj.noutwav)*ones((obj.nrow*obj.ncol),1);
            %%%%%%%%%%%%%%%%%%%

            
            
            %%%%%%%%%%%%%%%%%%%
            svarrshTV = reshape(svarTV,obj.nrow,obj.ncol,2);
            
            svarrshTV1 = svarrshTV(:,:,1);
            svarrshTV2 = svarrshTV(:,:,2);
            
            rvarrshTV1 = svarrshTV1(:,[1,1:end-1])+ svarrshTV1;
            rvarrshTV1(:,1) = svarrshTV1(:,1);
            rvarrshTV1(:,end) = svarrshTV1(:,1);
            
            rvarrshTV2 = svarrshTV2([1,1:end-1],:) + svarrshTV2;
            rvarrshTV2(1,:)= svarrshTV2(1,:);
            rvarrshTV2(end,:)= svarrshTV2(end-1,:);
            
            rvarrshTV = rvarrshTV1 + rvarrshTV2;
            rvarTV = rvarrshTV(:);
            %%%%%%%%%%%%%%%%%%%
            
            rvar = [rvarZ;rvarW;rvarTV];
%              rvar = rvarZ;
        end
                    
        % Get mask
        function [Isamp] = getIsamp(obj)
            Isamp = obj.Isamp;
        end  
       
    end
    
end

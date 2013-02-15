function test1d();

% parameters
N = 256; % signal length
% total_samp=100;% number of measurements
k = 30;% sparsity

% sparsifying transform
Psi = idct(eye(N));

% gaussian measurement matrix
Mtx1 = orth(randn(N,N));   

% partial Hadamard with random permuatation
Mtx2 = hadamard(N)/sqrt(N);
q2 = [1,randperm(N-1)+1];

p2 = randperm(N);

% block-diagonal Hadamard with random permutation;
blk_size = 8;
Mtx3=kron(eye(N/blk_size),hadamard(blk_size)/sqrt(blk_size));
q3 = randperm(N);
p3 = randperm(N);

% partial Hadamard with random sign reversal;
Mtx4 = hadamard(N)/sqrt(N);
p4 = 2*round(rand(N,1))-1;
Mtx4 = Mtx4*diag(p4);
q4 = randperm(N);


trial_num = 500;
success_Gauss = zeros(6,1);
success_Had_RP = zeros(6,1);
success_blkHad = zeros(6,1);
success_Had_RS = zeros(6,1);

for j =1:6
    samp_num = 60+(j-1)*10;
    % Gaussian Sampling Operator;
    A1 = Mtx1(1:samp_num,:);
    % Full WHT Sampling Operator with random permutation;
    A2 = Mtx2(q2(1:samp_num),p2);
     % Block WHT Sampling Operator with random permutation;
    A3 = Mtx3(q3(1:samp_num),p3);
    % Full WHT Sampling Operator with random sign flipping;
    A4 = Mtx4(q4(1:samp_num),:);
    
    message=sprintf('Sampling No.=%d',samp_num);
    disp(message);
    
    for i=1:trial_num

        % create a sparse signal in Psi domain
        alp = [randn(k,1); zeros(N-k,1)];

        p = randperm(N);
        alp = alp(p);
        x = Psi*alp;
        
        % observation
        y1 = A1*x;
        y2 = A2*x;
        y3 = A3*x;
        y4 = A4*x;
        
        % reconstruction using OMP alg.
        sigma = 0;
  
        alp1 = omp(x, A1*Psi, y1, k, sigma);
        xr1 = Psi*alp1;

        alp2 = omp(x, A2*Psi, y2, k, sigma);
        xr2 = Psi*alp2;

        alp3 = omp(x, A3*Psi, y3, k, sigma);
        xr3 = Psi*alp3;

        alp4 = omp(x, A4*Psi, y4, k, sigma);
        xr4 = Psi*alp4;

% if snr >50, it is considered as perfect reconstruction
        if snr(x,xr1)>50
            success_Gauss(j) = success_Gauss(j)+1;
        end
        
         if snr(x,xr2)>50
            success_Had_RP(j) = success_Had_RP(j)+1;
         end
        
          if snr(x,xr3)>50
            success_blkHad(j) = success_blkHad(j)+1;
          end
          
          if snr(x,xr4)>50
            success_Had_RS(j) = success_Had_RS(j)+1;
         end

    end
    
    success_Gauss(j) = success_Gauss(j)/trial_num;
    success_Had_RP(j) = success_Had_RP(j)/trial_num;
    success_blkHad(j) = success_blkHad(j)/trial_num;
    success_Had_RS(j) = success_Had_RS(j)/trial_num;

end 


figure;hold on;
plot([60:10:110], success_Gauss);
plot([60:10:110], success_Had_RP, 'k');
plot([60:10:110], success_blkHad, 'r');
plot([60:10:110], success_Had_RS, 'g');
legend('i.i.d Gaussian Operator','WHT256+Random Permutation','WHT8+Random Permutation','WHT256+Ramdom Sign Flipping');
xlabel('No. of Samples');
ylabel('Successful rate');
hold off;


%% this is the file used to test all 6 algorithms in the paper:
%  M. Yan, Y. Yang and S. Osher, Robust 1-bit compressive sensing using adaptive outlier pursuit. 
%  UCLA CAM report 11-71.
%%

clear all 
clc
%% Important parameters and functions
N = 1000;   % Signal dimension
M = 1000;   % Number of measurements
K = 10;     % Sparsity

L = 30;     % noise level


x0  = zeros(N,1);
rp  = randperm(N); 
x0(rp(1:K)) = randn(K,1);

A   = randn(M,N);
y0  = A * x0;  % Signal
b0  = sign(y0);

normalize = norm(x0);
x0  = x0/normalize;
y0  = y0/normalize;

rp      = randperm(M);
b       = b0;
b(rp(1:L)) = -b(rp(1:L));

%% Testing BIHT (6 versions)
x       = A'*b;
maxiter = 300;
alpha   = 1;

[x_1 iter_1] = BIHT(b, A, x, K, maxiter, alpha);

miter   = 1;

[x_3 iter_3 Loc_3] = BIHT_AOP(b, A, x, K, L, miter, maxiter, alpha);
[x_4 iter_4 Loc_4] = BIHT_AOP_flip(b, A, x, K, L, miter, maxiter, alpha);

x       = A'*b;
normalize = norm(x);
x       = x/normalize;  

alpha   = 1/M;

[x_2 iter_2] = BIHT2(b, A, x, K, maxiter, alpha);

[x_5 iter_5 Loc_5] = BIHT2_AOP(b, A, x, K, L, miter, maxiter, alpha);
[x_6 iter_6 Loc_6] = BIHT2_AOP_flip(b, A, x, K, L, miter, maxiter, alpha);



SNR     = 20*log10([1/norm(x_1/norm(x_1)-x0/norm(x0)), 1/norm(x_2/norm(x_2)-x0/norm(x0)), ...
                   1/norm(x_3/norm(x_3)-x0/norm(x0)),1/norm(x_4/norm(x_4)-x0/norm(x0)), ...
                   1/norm(x_5/norm(x_5)-x0/norm(x0)),1/norm(x_6/norm(x_6)-x0/norm(x0))])
HamErr  = [nnz(b0-sign(A*x_1)), nnz(b0-sign(A*x_2)), nnz(b0-sign(A*x_3)), ...
           nnz(b0-sign(A*x_4)), nnz(b0-sign(A*x_5)), nnz(b0-sign(A*x_6))]/M
AngErr  = [acos(x_1'*x0/norm(x_1)/norm(x0)), acos(x_2'*x0/norm(x_2)/norm(x0)), ...
           acos(x_3'*x0/norm(x_3)/norm(x0)), acos(x_4'*x0/norm(x_4)/norm(x0)), ...
           acos(x_5'*x0/norm(x_5)/norm(x0)), acos(x_6'*x0/norm(x_6)/norm(x0))]/pi       
HamDis  = [nnz(b-sign(A*x_1)), nnz(b-sign(A*x_2)), nnz(b-sign(A*x_3)), ...
           nnz(b-sign(A*x_4)), nnz(b-sign(A*x_5)), nnz(b-sign(A*x_6))]/M

scrsz = get(0,'ScreenSize');
figure('Position',[0.1*scrsz(3) 0.1*scrsz(4) 0.8*scrsz(3) 0.8*scrsz(4)])      
subplot(2,3,1)
plot(x0/norm(x0),'rx')
hold on
plot(x_1/norm(x_1),'bo')
legend('Original','BIHT','Location','NorthOutside')

subplot(2,3,2)
plot(x0/norm(x0),'rx')
hold on
plot(x_2/norm(x_2),'bo')
legend('Original','BIHT-L2','Location','NorthOutside')

subplot(2,3,3)
plot(x0/norm(x0),'rx')
hold on
plot(x_3/norm(x_3),'bo')
legend('Original','AOP','Location','NorthOutside')

subplot(2,3,4)
plot(x0/norm(x0),'rx')
hold on
plot(x_4/norm(x_4),'bo')
legend('Original','AOP-f','Location','NorthOutside')

subplot(2,3,5)
plot(x0/norm(x0),'rx')
hold on
plot(x_5/norm(x_5),'bo')
legend('Original','AOP-L2','Location','NorthOutside')

subplot(2,3,6)
plot(x0/norm(x0),'rx')
hold on
plot(x_6/norm(x_6),'bo')
legend('Original','AOP-L2-f','Location','NorthOutside')
% function [V, X0_V, img_V] = GF(X, Y, W, beta, lambda, eta, anneal, sigma0, iter_num, N0, is_grad, fgt, X0_X, img_X)
function [V] = GF(X, Y, W, beta, lambda, eta, anneal, sigma0, iter_num, N0, is_grad, fgt)
% Robust Point Set Registration Using Gaussian Mixture Models
% Zhao Ji, 2011-09-1
% Reference
% [Bin Jian 2010 PAMI]

n_ker = N0;
[N, D] = size(X); 
[M, D] = size(Y);
tmp_X = unique(X, 'rows'); idx = randperm(size(tmp_X,1)); idx = idx(1:min(n_ker,size(tmp_X,1))); ctrl_X=tmp_X(idx,:);
K=con_K(ctrl_X, ctrl_X, beta);
U = con_K(X, ctrl_X, beta);
% X0_U = con_K(X0_X, ctrl_X, beta);
% img_U = con_K(img_X, ctrl_X, beta);

x0 = zeros(n_ker*D, 1);%/N^2;
sigma = sigma0; %sum(sum((Y-X).^2))/(N*D) /36;

%%
options = optimset( 'display','off', 'MaxIter', 50);
% options = optimset( 'display','iter');
% options = optimset( 'display','iter','TolFun',1e-3,'TolX',1e-3);
if is_grad
    options = optimset(options, 'GradObj', 'on');
end
param = fminunc(@(x)costfun(x, X, Y, W, K, U, lambda, sigma, is_grad, fgt), x0, options);
for ii = 1:iter_num
    if ~fgt
        C = param(1:end);
        C = reshape(C, [n_ker D]);
        V = X + U*C;
        if eta
            W = comput_w(V, Y, eta);         
            W = W*N*M/sum(W(:));
        end
%         W = eye(size(X,1));
    end
    sigma = sigma*anneal;
    param = fminunc(@(x)costfun(x, X, Y, W, K, U, lambda, sigma, is_grad, fgt), param, options);
end

C = param(1:end);
C = reshape(C, [n_ker D]);
V = X + U*C;
% X0_V = X0_X + X0_U*C;
% img_V = img_X + img_U*C;
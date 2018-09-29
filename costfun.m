function [E, G] = costfun(param, X, Y, W, K, U, lambda, sigma, is_grad, fgt)

[N, D] = size(X);
[M, D] = size(Y);
N0 = size(U, 2);
C = reshape(param, [N0 D]);

E = lambda * trace(C'*K*C);

if fgt
    e = 10;      % Ratio of far field (default e = 10)
    K = sqrt(N); % Number of centers (default K = sqrt(Nx))
    p = 8;      % Order of truncation (default p = 8)

    [xc , A_k] = fgt_model((X+U*C)', ones(1,N), sigma, e,K,p);
    F = fgt_predict(Y' , xc , A_k , sigma, e);
    E = E-sum(F);
else
    T = X+U*C;
    F = pdist2(T,Y,'euclidean').^2;
    F =W.*exp(-F/sigma^2);
    E = E - sum(F(:));
end

%%
G = [];
if is_grad
    G = 2 * lambda * K * C;
    tmp = zeros(size(G));
    for m = 1:M
        tmp = tmp + bsxfun(@times, U, F(:,m))' * bsxfun(@minus, T, Y(m,:));
    end
    G = G + tmp*2/sigma^2;
    G = G(:);
end
% function [E, G] = costfun(param, X, Y, W, K, U, lambda, sigma, is_grad, fgt)
% 
% [N, D] = size(X);
% [M, D] = size(Y);
% N0 = size(U, 2);
% C = reshape(param, [N0 D]);
% 
% E = lambda * trace(C'*K*C);
% 
% if fgt
%     e = 10;      % Ratio of far field (default e = 10)
%     K = 20;sqrt(N); % Number of centers (default K = sqrt(Nx))
%     p = 8;      % Order of truncation (default p = 8)
% 
%     [xc, A_k] = fgt_model((X+U*C)', ones(1,N), sigma, e, K, p);
%     F = fgt_predict(Y', xc, A_k, sigma, e);
%     E = E-sum(F);
% else
%     for m=1:M
%         YM = repmat(Y(m,:),N,1);
%         WM = W(:,m);
%         V = YM-X-U*C;
%         F = WM .* exp(-sum(V.^2,2)/sigma^2);
%         E = E-sum(F);
%     end
% end
% 
% %%
% G = [];
% if is_grad
%     G = 2 * lambda * K * C;
%     for m=1:M
%         YM2 = repmat(Y(m,:),N,1);
%         WM2 = repmat(W(:,m),1,D);
%         V2 = -YM2+X+U*C;
%         A2 = repmat(exp(-sum(V2.^2,2)/sigma^2),1,D);
%         F2 = U'*(V2 .* WM2 .* A2)*2/sigma^2;
%         G = G+F2;
%     end
%     G = G(:);
% end
% 
% 

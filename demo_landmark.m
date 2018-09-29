
clear
% close all;
rand('state', 0);

if 1
    addpath('.\SC');
    addpath('.\FGT');
end

load('landmarks.mat');

j = 5;
X = landmarks(1:58,:,1);
Y = landmarks((j-1)*58+1:58*j,:,1);

[N, D] = size(X);
[M, D] = size(Y);

iter_num = 1;
anneal = 0.93;
sigma0 = 0.08;
N0 = 10;
eta = 1;
is_grad = 1;
fgt = 0;
beta = 0.5;
lambda = 0.01;

normalize = 1;
visualize = 1;

W = ones(size(X,1),size(Y,1));
if fgt
    is_grad = 0;
else
    % configuration of SC
    if eta
    W = comput_w(X, Y, eta); 
    W = W*N*M/sum(W(:));
    end
%     W = eye(size(X,1));
end

normal.xm=0; normal.ym=0;
normal.xscale=1; normal.yscale=1;
if normalize, [nX, nY, normal]=norm2(X,Y); end
% sigma0 = sqrt(sum(sum((nY(1:N,:)-nX).^2))/(N*D))/3;

tic;
V = GF(nX, nY, W, beta, lambda, eta, anneal, sigma0, iter_num, N0, is_grad, fgt);
t=toc

if normalize, V=V*normal.yscale+repmat(normal.ym,size(X,1),1); end 

% mea = sum(sqrt(sum((y2(1:size(V,1),:)-V).^2,2)))/size(V,1)
   
if visualize
  figure;
  plot(X(:,1),X(:,2),'b+',Y(:,1),Y(:,2),'ro')
  title('before');
  hold off
  drawnow	
  figure;
  plot(V(:,1),V(:,2),'b+',Y(:,1),Y(:,2),'ro')
  title('after');
  hold off
  drawnow	
end
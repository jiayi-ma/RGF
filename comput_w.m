function W = comput_w(X, Y, eta)

mean_dist_global=[]; % use [] to estimate scale from the data
nbins_theta=12;
nbins_r=5;
nsamp1=size(X,1);
nsamp2=size(Y,1);
r_inner=1/8;
r_outer=2;

% if visualize
%    figure(21)
%    plot(X(:,1),X(:,2),'b+',Y(:,1),Y(:,2),'ro')
%    title(['original pointsets (nsamp1=' int2str(nsamp1) ', nsamp2=' ...
%        int2str(nsamp2) ')'])
%    drawnow
% end

out_vec_1=zeros(1,nsamp1); 
out_vec_2=zeros(1,nsamp2);

[BH1,mean_dist_1]=sc_compute(X',zeros(1,nsamp1),mean_dist_global,nbins_theta,nbins_r,r_inner,r_outer,out_vec_1);
[BH2,mean_dist_2]=sc_compute(Y',zeros(1,nsamp2),mean_dist_1,nbins_theta,nbins_r,r_inner,r_outer,out_vec_2);
costmat=hist_cost_2(BH1,BH2);

W = exp(-eta*costmat);
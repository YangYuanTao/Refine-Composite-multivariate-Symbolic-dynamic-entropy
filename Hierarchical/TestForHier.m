clear;clc;close;

N=1024;
Data=randn(1,N);
Normal_Data=randn(1,N);

% The hierarchical analysis below is modified by Yang Yuantao and Li Yongbo .
% hierarchical sample entropy
m=2;r=0.15;t=1;n=3;
HSE=HierSE(Data,m,r,t,n);
% hierarchical fuzzy entropy
m=2;r=0.15;n_FE=2;t=1;n=3;
HFE=HierFE(Data,m,r,n_FE,t,n);
% hierarchical dispersion entropy
m=2;nc=6;t=1;n=3;
HDE=HierDE(Data,m,nc,t,n);
% hierarchical permutation entropy
m=3;t=1;n=3;
HPE=HierPE_mu(Data,m,t,n);
% hierarchical symbolic dynamic entropy
numSymbol=6;Depth=3;n=3;
HSDE=HierSymbolicDynamicEntropy(Normal_Data,Data,numSymbol,Depth,n);


figure(1)
x=1:8;
hold on
plot(x,HSE,'.-','Displayname','HSE','markersize',12)
plot(x,HFE,'.-','Displayname','HFE','markersize',12)
plot(x,HDE,'.-','Displayname','HDE','markersize',12)
plot(x,HPE,'.-','Displayname','HPE','markersize',12)
plot(x,HSDE,'.-','Displayname','HSDE','markersize',12)
hold off
legend('show','Location','Best')
xlim([0,9])
xlabel('Hier node')
ylabel('Entropy value')
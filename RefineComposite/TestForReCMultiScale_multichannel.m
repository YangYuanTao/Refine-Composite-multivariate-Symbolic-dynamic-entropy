clear;clc;close;

N=1024;
Data=randn(3,N);
Normal_Data=randn(3,N);

% hierarchical sample entropy
m=2;r=0.15;t=1;Scale=8;
RCmvMSE=RCmvMSE_mu(Data,m,r,t,Scale);
% hierarchical fuzzy entropy
m=2;r=0.15;n_FE=2;t=1;Scale=8;
RCmvMFE=RCmvMFE_mu(Data,m,r,n_FE,t,Scale);
% hierarchical permutation entropy
m=3;t=1;Scale=8;
RCmvMPE=RCmvMPE_mu(Data,m,Scale);

% hierarchical symbolic dynamic entropy
numSymbol=6;Depth=3;Scale=8;
RCmvMSDE=RCmvMSDE_MEP(Normal_Data,Data,numSymbol,Depth,Scale);


figure(1)
x=1:8;
hold on
plot(x,RCmvMSE,'.-','Displayname','RCmvMSE','markersize',12)
plot(x,RCmvMFE,'.-','Displayname','RCmvMFE','markersize',12)
plot(x,RCmvMPE,'.-','Displayname','RCmvMPE','markersize',12)
plot(x,RCmvMSDE,'.-','Displayname','RCmvMSDE','markersize',12)
hold off
legend('show','Location','Best')
xlim([0,9])
xlabel('Scale')
ylabel('Entropy value')
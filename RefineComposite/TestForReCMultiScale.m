clear;clc;close;

N=1024;
Data=randn(1,N);
Normal_Data=randn(1,N);

% hierarchical sample entropy
m=2;r=0.15;t=1;Scale=8;
RCMSE=RCmvMSE_mu(Data,m,r,t,Scale);
% hierarchical fuzzy entropy
m=2;r=0.15;n_FE=2;t=1;Scale=8;
RCMFE=RCmvMFE_mu(Data,m,r,n_FE,t,Scale);
% hierarchical dispersion entropy
m=2;nc=6;t=1;Scale=8;
RCMDE=RCMDE(Data,m,nc,t,Scale);
% hierarchical permutation entropy
m=3;t=1;Scale=8;
RCMPE=RCmvMPE_mu(Data,m,Scale);

% hierarchical symbolic dynamic entropy
numSymbol=6;Depth=3;Scale=8;
RCMSDE=RCmvMSDE_MEP(Normal_Data,Data,numSymbol,Depth,Scale);


figure(1)
x=1:8;
hold on
plot(x,RCMSE,'.-','Displayname','RCMSE','markersize',12)
plot(x,RCMFE,'.-','Displayname','RCMFE','markersize',12)
plot(x,RCMDE,'.-','Displayname','RCMDE','markersize',12)
plot(x,RCMPE,'.-','Displayname','RCMPE','markersize',12)
plot(x,RCMSDE,'.-','Displayname','RCMSDE','markersize',12)
hold off
legend('show','Location','Best')
xlim([0,9])
xlabel('Scale')
ylabel('Entropy value')
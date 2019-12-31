clear;clc;close;

N=1024;
Data=randn(1,N);
Normal_Data=randn(1,N);

% hierarchical sample entropy
m=2;r=0.15;t=1;Scale=8;
MMSE=Modified_mvMSE_mu(Data,m,r,t,Scale);
% hierarchical fuzzy entropy
m=2;r=0.15;n_FE=2;t=1;Scale=8;
MMFE=Modified_mvMFE_mu(Data,m,r,n_FE,t,Scale);
% hierarchical dispersion entropy
m=2;nc=6;t=1;Scale=8;
MMDE=Modified_MDE(Data,m,nc,t,Scale);
% hierarchical permutation entropy
m=3;t=1;Scale=8;
MMPE=Modified_mvMPE_mu(Data,m,Scale);

% hierarchical symbolic dynamic entropy
numSymbol=6;Depth=3;Scale=8;
MMSDE=Modifiied_MSymbolicDynamicEntropy(Normal_Data,Data,numSymbol,Depth,Scale);


figure(1)
x=1:8;
hold on
plot(x,MMSE,'.-','Displayname','MMSE','markersize',12)
plot(x,MMFE,'.-','Displayname','MMFE','markersize',12)
plot(x,MMDE,'.-','Displayname','MMDE','markersize',12)
plot(x,MMPE,'.-','Displayname','MMPE','markersize',12)
plot(x,MMSDE,'.-','Displayname','MMSDE','markersize',12)
hold off
legend('show','Location','Best')
xlim([0,9])
xlabel('Scale')
ylabel('Entropy value')
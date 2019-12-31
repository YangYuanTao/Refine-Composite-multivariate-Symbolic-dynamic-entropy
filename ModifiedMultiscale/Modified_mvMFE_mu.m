function Out_mvMFE=Modified_mvMFE_mu(X,sm,sr,n,stau,Scale)
% This function calculates multivariate multiscale fuzzy entropy (mvMFE) whose
% coarse-graining is based on mean
%
%
% Inputs:
% X: multivariate signal - a matrix of size nvar (the number of channels) x nsamp (the number of sample points for each channel)
% sm: scalar embedding value
% sr: scalar threshold value (it is usually equal to 0.15)
% n: fuzzy power (it is usually equal to 2)
% stau: scalar time lag  value (it is usually equal to 1)
% Scale: the number of scale factors
%
% Output:
% Out_mvMFE: a scalar quantity
%
% Ref:
% [1] H. Azami and J. Escudero, "Refined Composite Multivariate Generalized Multiscale Fuzzy Entropy:
% A Tool for Complexity Analysis of Multichannel Signals", Physica A, 2016.
%
% If you use the code, please make sure that you cite reference [1].
%
% Hamed Azami and Javier Escudero Rodriguez
% hamed.azami@ed.ac.uk and javier.escudero@ed.ac.uk
%
%  10-June-16


%%
% Because multi-channel signals may have different amplitude ranges, the distances calculated on embedded vectors may be biased
%toward the largest amplitude ranges variates. For this reason, we scale all the data channels to the same amplitude range and
% we normalize each data channel to unit standard deviation so that the total variation becomes equal to the number of channels
% or variables [1].

X=zscore(X');
r = sr*(sum(std(X)));
M = sm*ones(1,size(X,2)); % M: embedding vector
tau=ones(1,size(X,2))*stau; % tau: time lag vector
X=X';

Out_mvMFE(1)=mvFE(X,M,r,n,tau);
for k=2:Scale
    Xs = Multi(X,k);
    Out_mvMFE(k) = mvFE(Xs,M,r,n,tau);
end

end


function M_Data = Multi(Data,S)
%  generate the consecutive coarse-grained time series based on mean
%  Input:   Data: time series;
%           S: the scale factor

% Output:
%           M_Data: the coarse-grained time series at the scale factor S
L = size(Data,2);
M_Data=NaN*ones(size(Data,1),L-S+1);
for j=1:size(Data,1)
    for i=1:L-S+1
        M_Data(j,i) = mean(Data(j,(i:i+S-1)));
    end
end
end

function [mvFE_Out,fi_m1,fi_m2]=mvFE(X,M,r,n,tau)
%
% This function calculates multivariate fuzzy entropy (mvFE) of a multivariate signal
%
% Inputs:
%
% X: multivariate signal - a matrix of size nvar (the number of channels) x nsamp (the number of sample points for each channel)
% M: embedding vector
% r: scalar threshold
% n: fuzzy power (it is usually equal to 2)
% tau: time lag vector
%
% Outputs:
%
% mvFE_Out: scalar quantity - the mvFE of X
% fi_m1: scalar quantity - the global quantity in dimension m
% fi_m2 : scalar quantity - the global quantity in dimension m+1
%
%
% Ref:
% [1] H. Azami and J. Escudero, "Refined Composite Multivariate Generalized Multiscale Fuzzy Entropy:
% A Tool for Complexity Analysis of Multichannel Signals", Physica A, 2016.
%
%
% If you use the code, please make sure that you cite reference [1].
%
% Hamed Azami and Javier Escudero Rodriguez
% hamed.azami@ed.ac.uk and javier.escudero@ed.ac.uk
%
%  10-June-16


max_M=max(M);
max_tau=max(tau);
nn=max_M*max_tau;


[nvar,nsamp]=size(X);
N=nsamp-nn;
A=embd(M,tau,X); % all the embedded vectors are created based on the Takens embedding Theorem - the code is publicly-available in
% "http://www.commsp.ee.ic.ac.uk/~mandic/research/Complexity_Stuff.htm" prepared by Prof. Mandic's group

y=pdist(A,'chebychev'); %infinite norm is calculated between all possible pairs

% For a given threshold r and fuzzy power n, define a global quantity fi_m1
y=exp((-y.^n)/r);
fi_m1=sum(y)*2/(N*(N-1));

clear yn y v1 A;


%% extend the dimensionality of the multivariate delay vector from m to
% m+1 in nvar differnet ways
M=repmat(M,nvar,1);
I=eye(nvar);
M=M+I;

B=[];

for h=1:nvar
    Btemp=embd(M(h,:),tau,X);
    B=vertcat(B,Btemp);% all the delay embedded vectors of all the subspaces of dimension m+1 is concatenated into a single matrix
    Btemp=[];
end

z=pdist(B,'chebychev'); %now comparison is done between subspaces
z=exp((-z.^n)/r);
fi_m2=sum(z)*2/(nvar*N*(nvar*N-1));

%% calculating MVFE

mvFE_Out=log(fi_m1/fi_m2);

end


function A=embd(M,tau,ts)
% This function creates multivariate delay embedded vectors with embedding
% vector parameter M and time lag vector parameter tau.
% M is a row vector [m1 m2 ...mnvar] and tau is also a row vector [tau1 tau2....taunvar] where nvar is the 
% number of channels; 
% ts is the multivariate time series-a matrix of size nvarxnsamp;

[nvar,nsamp]=size(ts);
A=[];

for j=1:nvar
for i=1:nsamp-max(M)*max(tau)
temp1(i,:)=ts(j,i:tau(j):i+(M(j)-1)*tau(j));
end
A=horzcat(A,temp1);
temp1=[];
end
end
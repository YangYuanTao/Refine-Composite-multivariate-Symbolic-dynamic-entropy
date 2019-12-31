function HFE=HierFE(data,m,r,n_FE,t,n)
%  Calculate the Hierarchical Fuzzy Entropy (HFE)
%  Input:   data: time series;
%           m: embedding dimension;
%           r: std of data;
%           t: time delay
%           n: floor of the Hierarchical Entropy;
%  Output: 
%           HFE: hierarchical fuzzy entropy.
%  code is arranged by yyt in 2019.01     yangyuantaohit@163.com
if size(data,1)==1
   data=data(:);
end
HFE =zeros(1,2^n);
for i=0:(2^n-1)
    datani = HierarchicalEn(n,i,data);

    % datani= n * 1;
    if size(datani,2)==1
        datani=datani';
    end
    sr=r*std(datani);
    FE = mvFE(datani,m,sr,n_FE,t);
    HFE(i+1) =FE;
end
end

function datani = HierarchicalEn(n,i,data)

imatrix=iarray(n,i);
for j=1:n
   [Q0,Q1]=QOperator(j,data);
    if imatrix(j)==1
        data = Q1*data;
    else
        data = Q0*data;
    end
end

datani=data;
end

function [Q0,Q1]=QOperator(n,data)
N=length(data);

Q0=zeros(N-2^(n-1),N);
for i=1:N-2^(n-1)
    Q0(i,i)=1/2;
    Q0(i,i+2^(n-1))=1/2;
end

Q1=zeros(N-2^(n-1),N);
for i=1:N-2^(n-1)
    Q1(i,i)=1/2;
    Q1(i,i+2^(n-1))=-1/2;
end

end

function imatrix=iarray(n,i)
istr=dec2bin(i);
if length(istr)==n
    imatrix=zeros(1,n);
    for j=1:n
        imatrix(j)=str2double(istr(j));
    end    
else    
    hh=n-length(istr);
    imatrix=zeros(1,n);
    imatrix(1:hh)=0;
    for j=(hh+1):n
    imatrix(j)=str2double(istr(j-hh));
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
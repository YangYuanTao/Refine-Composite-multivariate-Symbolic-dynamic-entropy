function HSE=HierSE(data,m,r,t,n)
%  Calculate the Hierarchical Sample Entropy (HSE)
%  Input:   data: time series;
%           m: embedding dimension;
%           r: std of data;
%           t: time delay;
%           n: floor of the Hierarchical Entropy;
%  Output: 
%           HSE: hierarchical sample entropy.
%  code is arranged by yyt in 2019.01     yangyuantaohit@163.com
if size(data,1)==1
   data=data(:);
end
HSE =zeros(1,2^n);
for i=0:(2^n-1)
    datani = HierarchicalEn(n,i,data);
    if size(datani,2)==1
        datani=datani';
    end
    SE = mvSE(datani,m,r,t);
    HSE(i+1) =SE;
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

function [Out_mvSE,p1,p2]=mvSE(X,M,r,tau)
%
% This function calculates multivariate sample entropy (mvSE) of a multivariate signal
%
% Inputs:
%
% X: multivariate signal - a matrix of size nvar (the number of channels) x nsamp (the number of sample points for each channel)
% M: embedding vector
% r: scalar threshold
% tau: time lag vector
%
% Outputs:
%
% e: scalar quantity - the mvSE of X
% p1: scalar quantity - the probability that any two composite delay vectors are similar in dimension m
% p2 : scalar quantity - the probability that any two composite delay vectors are similar in dimension m+1
%
% The code is based on the publicly-available code in
% "http://www.commsp.ee.ic.ac.uk/~mandic/research/Complexity_Stuff.htm"
% prepared by Prof. Mandic's group
%
% Ref:
% [1] M. U. Ahmed and D. P. Mandic, "Multivariate multiscale entropy
% analysis", IEEE Signal Processing Letters, vol. 19, no. 2, pp.91-94.2012
% [2] H. Azami and J. Escudero, "Refined Composite Multivariate Generalized Multiscale Fuzzy Entropy:
% A Tool for Complexity Analysis of Multichannel Signals", Physica A, 2016.
%
%
% If you use the code, please make sure that you cite references [1] and [2].
%
% Hamed Azami and Javier Escudero Rodriguez
% hamed.azami@ed.ac.uk and javier.escudero@ed.ac.uk
%
%  10-June-16
r=r*std(X);
mm=max(M);
mtau=max(tau);
nn=mm*mtau;


[nvar,nsamp]=size(X);
N=nsamp-nn;
A=embd(M,tau,X);%all the embedded vectors are created
y=pdist(A,'chebychev');%infinite norm is calculated between all possible pairs
[r1,c1,v1]= find(y<=r);% threshold is implemented
p1=numel(v1)*2/(N*(N-1));%the probability that two templates of length m are closed within the tolerance r
clear  y r1 c1 v1 A;

M=repmat(M,nvar,1);
I=eye(nvar);
M=M+I;

B=[];

% number of match templates of length m+1 closed within the tolerance r where m=sum(M) is calculated afterwards
for h=1:nvar
    Btemp=embd(M(h,:),tau,X);
    B=vertcat(B,Btemp);% all the delay embedded vectors of all the subspaces of dimension m+1 is concatenated into a single matrix
    Btemp=[];
end
z=pdist(B,'chebychev'); %now comparison is done between subspaces
[r2,c2,v2]= find(z<=r);
p2=numel(v2)*2/(nvar*N*(nvar*N-1));
clear  z r2 c2 v2 B;

Out_mvSE=log(p1/p2);
end

function A=embd(M,tau,ts)
% This function creates multivariate delay embedded vectors with embedding
% vector parameter M and time lag vector parameter tau.
% M is a row vector [m1 m2 ...mnvar] and tau is also a row vector [tau1 tau2....taunvar] where nvar is the 
% number of channels; 
% ts is the multivariate time series-a matrix of size nvarxnsamp;

% Ref: M. U. Ahmed and D. P. Mandic, "Multivariate multiscale entropy
% analysis", IEEE Signal Processing Letters, vol. 19, no. 2, pp.91-94.2012

[nvar,nsamp]=size(ts);
A=[];

for j=1:nvar
for i=1:nsamp-max(M)
    temp1(i,:)=ts(j,i:tau(j):i+M(j)-1);
end
A=horzcat(A,temp1);
temp1=[];
end
end
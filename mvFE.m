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

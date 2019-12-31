
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
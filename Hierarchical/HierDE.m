function HDE=HierDE(Data,m,nc,t,n)
%  Calculate the Hierarchical Dispersion Entropy
%  Input:   data: time series;
%           m: embedding dimension;
%           nc: number of classes;
%           t: time delay ;
%           n: the scale of Hierarchical;
%  Output: 
%           HDE: the Hierarchical Dispersion Entropy.
%  code is arranged by yyt in 2019.01     yangyuantaohit@163.com
if size(Data,1)==1
   Data=Data(:);
end
HDE=zeros(1,2^n);
for i=0:(2^n-1)
    Hier_Data = HierEn(n,i,Data);
    if size(Hier_Data,2)==1
        Hier_Data=Hier_Data';
    end
    HDE(i+1)= DisEn_NCDF(Hier_Data,m,nc,t);
end
end

function [datani] = HierEn(n,i,data)
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

function [Out_DisEn, npdf]=DisEn_NCDF(x,m,nc,tau)
%
% This function calculates dispersion entropy (DisEn) of a univariate
% signal x, using normal cumulative distribution function (NCDF)
%
% Inputs:
%
% x: univariate signal - a vector of size 1 x N (the number of sample points)
% m: embedding dimension
% nc: number of classes (it is usually equal to a number between 3 and 9 - we used c=6 in our studies)
% tau: time lag (it is usually equal to 1)
%
% Outputs:
%
% Out_DisEn: scalar quantity - the DisEn of x
% npdf: a vector of length nc^m, showing the normalized number of disersion patterns of x
%
% Ref:
%
% [1] H. Azami, M. Rostaghi, D. Abasolo, and J. Escudero, "Refined Composite Multiscale Dispersion Entropy and its Application to Biomedical
% Signals", IEEE Transactions on Biomedical Engineering, 2017.
% [2] M. Rostaghi and H. Azami, "Dispersion Entropy: A Measure for Time-Series Analysis", IEEE Signal Processing Letters. vol. 23, n. 5, pp. 610-614, 2016.
%
% If you use the code, please make sure that you cite references [1] and [2].
%
% Hamed Azami, Mostafa Rostaghi, and Javier Escudero Rodriguez
% hamed.azami@ed.ac.uk, rostaghi@yahoo.com, and javier.escudero@ed.ac.uk
%
%  20-January-2017
%%


N=length(x);

sigma_x=std(x);
mu_x=mean(x);

y=normcdf(x,mu_x,sigma_x);

for i_N=1:N
    if y(i_N)==1
        y(i_N)=1-(1e-10);
    end
    
     if y(i_N)==0
        y(i_N)=(1e-10);
     end
end

z=round(y*nc+0.5);

all_patterns=[1:nc]';

for f=2:m
    temp=all_patterns;
    all_patterns=[];
    j=1;
    for w=1:nc
        [a,b]=size(temp);
        all_patterns(j:j+a-1,:)=[temp,w*ones(a,1)];
        j=j+a;
    end
end

for i=1:nc^m
    key(i)=0;
    for ii=1:m
        key(i)=key(i)*10+all_patterns(i,ii);
    end
end


embd2=zeros(N-(m-1)*tau,1);
for i = 1:m, 
    embd2=[z(1+(i-1)*tau:N-(m-i)*tau)]'*10^(m-i)+embd2;
end


pdf=zeros(1,nc^m);

for id=1:nc^m
    [R,C]=find(embd2==key(id));
    pdf(id)=length(R);
end

npdf=pdf/(N-(m-1)*tau);
p=npdf(npdf~=0);
Out_DisEn = -sum(p .* log(p));
end
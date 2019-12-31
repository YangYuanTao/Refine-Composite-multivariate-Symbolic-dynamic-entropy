function HDE=HierDE(Data,m,nc,t,n)
%  Calculate the Hierarchical Dispersion Entropy
%  Input:   data: time series;
%           m: embedding dimension;
%           t: time delay of permuation entropy;
%           nc: number of classes;
%           n: the scale of Hierarchical;
%  Output: 
%           HSDE: the Hierarchical Symbolic Dynamic Entropy.
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
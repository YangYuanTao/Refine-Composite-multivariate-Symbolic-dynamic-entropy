function HSE=HierSE(data,m,r,n)
%  Calculate the Hierarchical Sample Entropy (HSE)
%  Input:   data: time series;
%           m: embedding dimension;
%           r: std of data;
%           n: floor of the Hierarchical Entropy;
%  Output: 
%           HSE: hierarchical sample entropy.
%  code is arranged by yyt in 2016.07     yangyuantaohit@163.com
if size(data,1)==1
   data=data(:);
end
HSE =zeros(1,2^n);
for i=0:(2^n-1)
    datani = HierarchicalEn(n,i,data);
%     SE = SampleEntropy(datani,m,r);
    if size(datani,2)==1
        datani=datani';
    end
    SE = mvSE(datani,m,r,1);
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
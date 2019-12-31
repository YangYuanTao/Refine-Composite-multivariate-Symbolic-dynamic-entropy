function MHMorphMatrix=ModifiedHierMorpyMartix_NCDF(data,Depth,numSymbol,n)
% Compute the Probability Vector
% using normal cumulative distribution function (NCDF)
%
% Input:   
%         Data:        the raw data - a vector of size 1 x N (the number of sample points)
%         Depth:       the embedding demension
%         numSymbol:   the number of symbol
%           n:         floor of the Hierarchical analysis;
%  Output: 
%         MHMorphMatrix: Modified hierarchical MorphMatrix.
%  code is arranged by yyt in 2018.05     yangyuantaohit@163.com

% time series data is a N*1 vector
[a,~]=size(data);
if a==1
   data=data';
end 

MHMorphMatrix =zeros(numSymbol^(Depth+1),2^n);

for i=0:(2^n-1)
    datani = HierarchicalData(n,i,data);
    MHMorphMatrix(:,i+1) = MorphMatrix_NCDF(datani',Depth,numSymbol);
end
end

function datani = HierarchicalData(n,i,data)

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

function MorphMatrix = MorphMatrix_NCDF(Data,Depth,numSymbol)
% Compute the Probability Vector
% using normal cumulative distribution function (NCDF)
%
% Input:   
%         Data:        the raw data - a vector of size 1 x N (the number of sample points)
%         Depth:       the embedding demension
%         numSymbol:   the number of symbol
% Output:   
%         PVec:        the Probability Vector of data
%  code is arranged by yyt in 2017.07   yangyuantaohit@163.com

if min(size(Data))~=1
    Data=Data(:);
end

L = length(Data);
y = normcdf(Data,mean(Data),std(Data)); 

epsilon=1e-10;
for i=1:L
    if y(i)>1-epsilon
        y(i)=1-epsilon;
    end
    if y(i)<epsilon
        y(i)=epsilon;
    end
end

symbolSeq=round(y*numSymbol+0.5);

%% generate the embedding vector and compute the pVector 
numStates=numSymbol^Depth;
MorphMatrix = zeros(numStates,numSymbol);
% PVec = zeros(1,numStates);
for jj = 1:L-Depth
    WindowSet=symbolSeq(jj:jj+Depth);
    RowNum=1;
    for i=1:Depth
        RowNum=RowNum+(WindowSet(i)-1)*numSymbol^(Depth-i);
    end
    ColumnNum=WindowSet(1+Depth);
%     PVec(RowNum)=1+PVec(RowNum);
    MorphMatrix(RowNum,ColumnNum)=1+MorphMatrix(RowNum,ColumnNum);   
end
% PVec = PVec/sum(PVec);  % normalizing the computed vector

% normalize each row of Matrix to make it a stochastic matrix
for i=1:numStates
    rowSum = sum(MorphMatrix(i,:));
    if(rowSum~=0)
        MorphMatrix(i,:) = MorphMatrix(i,:)/rowSum;     
    end
end
MorphMatrix=MorphMatrix';
MorphMatrix=MorphMatrix(:);
end
function [HSDE]=HierarchicalSymbolicDynamicEntropy(PartitionData,Data,numSymbol,Depth,n)
%  compute the Hierarchical Symbolic Dynamic Entropy
%  Input:   PartitionData: the 'normal' data, either a vector or a matrix; 
%           (to generate the partition vector by Maximun Entropy  Partition)
%           Data: the time series, a vector;
%           numSymbol:  the number of the symbol kinds
%           Depth: the depth of the state;
%           n: the scale of Hierarchical;
%  Output: 
%           HSDE: the Hierarchical Symbolic Dynamic Entropy.
HSDE=zeros(1,2^n);
for i=0:(2^n-1)
    Hier_PartitionData = HierEn(n,i,PartitionData);
    partition = PartitionGeneration(Hier_PartitionData,numSymbol);
    
    Hier_Data = HierEn(n,i,Data);
    symbolSeq=SymbolGeneration(Hier_Data,partition);
    HSDE(i+1)=NormalSymbolicDynamicEntropy(symbolSeq, numSymbol, Depth);
end

end


function partition = PartitionGeneration(PartitionData,numSymbol)
%  generate the partition vector
%  Input:    PartitionData: time series¡¡as a column Vector
%            numSymbol: the number of the symbol kinds
%  Output:   
%            partition£ºthe patition value of the time series
[m,n] = size(PartitionData);
data = reshape(PartitionData,m*n,1);  %change into long vector
partition = maxEntropyPartition(data,numSymbol);
end

function partition = maxEntropyPartition(data,numSymbol)
%  generate the partition vector
%  Input:    PartitionData: time series¡¡as a column Vector
%            numSymbol: the number of the symbol kinds
%  Output:   
%            partition£ºthe patition value of the time series
x = sort(data);
k = max(size(x));
partition = zeros(1, numSymbol-1);
for i=1:numSymbol-1
    partition(1,i) = x(floor(i*k/numSymbol));
end
end

function symbolSeq=SymbolGeneration(Data,partition)
%  Calculate the symbol of the Data 
%  Input:    Data: time series¡¡as a column Vector
%            numSymbol: the number of the symbol kinds
%  Output:   
%            symbol£ºthe symbol sequence of the time series
[a,~]=size(Data);
if a==1
   Data=Data';
end 

partition=[partition inf];
symbolSeq=zeros(length(Data),1);
for ii = 1:length(Data)
    for s=1:length(partition)
        if partition(s)>Data(ii)
            symbolSeq(ii)=s;
            break
        end        
    end
end
end

function  [apen]=NormalSymbolicDynamicEntropy(symbolSeq, numSymbol, Depth)
%  Calculate the Symbol Entropy of the Data 
%  Input:    symbolSeq: time series¡¡as a symbol sequence;
%            numSymbol: the number of the symbol kinds;
%            Depth: the depth of the state;
%  Output:   
%            apen£ºthe symbol entropy of the symbol sequence.
% Note that:
% numStates=numSymbol^Depth;
% stateProbVecStore = zeros(numStates,1);
% MorphMatrix = zeros(numStates,numSymbol);
[MorphMatrix,stateProbVecStore] = AnalyzeSymbolSeq(symbolSeq, numSymbol, Depth);

apen=0;
A=[];
for i=1:length(stateProbVecStore) 
    for j=1:numSymbol
       if MorphMatrix(i,j)~=0
           A=[A,stateProbVecStore(i)*MorphMatrix(i,j)];
       end
    end
end 
apen=apen-sum(A.*log(A));
apen=apen/log(numSymbol^(Depth+1)); % normalizing the computed entropy
end

function [MorphMatrix, PVec] = AnalyzeSymbolSeq(symbolSeq,numSymbol,Depth)
numStates=numSymbol^Depth;
MorphMatrix = zeros(numStates,numSymbol);
PVec = zeros(1,numStates);

N=length(symbolSeq);
for jj = 1:N-Depth
    WindowSet=symbolSeq(jj:jj+Depth);
    RowNum=1;
    for i=1:Depth
        RowNum=RowNum+(WindowSet(i)-1)*numSymbol^(Depth-i);
    end
    ColumnNum=WindowSet(1+Depth);
    PVec(RowNum)=1+PVec(RowNum);
    MorphMatrix(RowNum,ColumnNum)=1+MorphMatrix(RowNum,ColumnNum);   
end
PVec = PVec/sum(PVec);  % normalizing the computed vector

%normalize each row of Matrix to make it a stochastic matrix
for i=1:numStates
    rowSum = sum(MorphMatrix(i,:));
    if(rowSum~=0)
        MorphMatrix(i,:) = MorphMatrix(i,:)/rowSum;     
    end
end
end

function [datani] = HierEn(n,i,data)
%%%%Hierarchical entropy Ñù±¾ÉÌ
% N=length(data);
% l=nextpow2(N);
% datani=zeros(1,N/(2^n));
imatrix=iarray(n,i);
for j=1:n
   [Q0,Q1]=QOperator(data);
    if imatrix(j)==1
        data = Q1*data;
    else
        data = Q0*data;
    end
end
datani=data;
end

function [Q0,Q1]=QOperator(data)
N=length(data);
l=nextpow2(N);

Q0=zeros(2^(l-1),2^l);
for i=1:2^(l-1)
    Q0(i,2*i-1)=1/2;
    Q0(i,2*i)=1/2;
end

Q1=zeros(2^(l-1),2^l);
for i=1:2^(l-1)
    Q1(i,2*i-1)=1/2;
    Q1(i,2*i)=-1/2;
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
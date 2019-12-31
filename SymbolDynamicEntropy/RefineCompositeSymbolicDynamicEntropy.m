function [RCSDE]=RefineCompositeSymbolicDynamicEntropy(PartitionData,Data,numSymbol,Depth,Scale)
%  compute the Refine Composite Symbolic Dynamic Entropy
%  Input:   PartitionData: the 'normal' data, either a vector or a matrix; 
%           (to generate the partition vector by Maximun Entropy  Partition)
%           Data: the time series, a vector;
%           numSymbol:  the number of the symbol kinds
%           Depth: the depth of the state;
%           Scale: the scale factor;
%  Output: 
%           RCSDE: the Refine Composite Symbolic Dynamic Entropy.
RCSDE=zeros(Scale,1);
for i=1:Scale
    Multi_PartitionData = Multi(PartitionData,Scale,2);
    partition = PartitionGeneration(Multi_PartitionData,numSymbol);
    MSDE=zeros(1,i);
    for j=1:i
        Multi_Data=Multi(Data(i:end),i,2);
        symbolSeq=SymbolGeneration(Multi_Data,partition);
        MSDE(j)=NormalSymbolicDynamicEntropy(symbolSeq, numSymbol, Depth);
    end
    RCSDE(i)=mean(MSDE);
end
end

function M_Data = Multi(data,Scale,flag)
%  generate the consecutive coarse-grained time series
%  Input:   data: time series;
%           scale: the scale factor;
%  Output: 
%           M_Data: the coarse-grained time series at the scale factor S.
[m,n]=size(data);
for k=1:n
    %% Multi-    
    if flag==1
        J = fix(m/Scale);
        M_Data=zeros(J,n);
        for i=1:J
            M_Data(i,k) = mean(data((i-1)*Scale+1:i*Scale,k));
        end
    elseif flag==2
    %% Modified Multi-
        J = m-Scale+1;
         M_Data=zeros(J,n);
         for i=1:J
            M_Data(i,k) = mean(data(i:i + Scale - 1,k));
         end
    end   
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
function  [apen]=SymbolEntropy(symbolSeq, numSymbol, Depth)
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
%
%¡¡code is arranged by yyt in 2016.06     yangyuantaohit@163.com
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
apen=apen/log(numSymbol^(Depth+1));
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
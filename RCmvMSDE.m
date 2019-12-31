function Out_mvMSDE=RCmvMSDE(Data,numSymbol,Depth,scale)
%  Inputs:   
%           Data: the time series, a matrix;
%           numSymbol:  the number of the symbol kinds
%           Depth: the depth of the state;
%           Scale: the scale factor;
%  Output: 
%           Out_mvMSDE: the Refine Composite SpatioTemporal Entropy
%  code is arranged by yyt in 2018.08   yangyuantaohit@163.com
Out_mvMSDE=NaN*ones(1,scale);
for i=1:scale
    temp_MorphMatrix=NaN*ones(numSymbol^Depth,numSymbol,i);
    temp_Pvector=NaN*ones(numSymbol^Depth,i);
    for j=1:i
        Multi_Data=Multi(Data(:,j:end),i);
        [~,temp_MorphMatrix(:,:,j),temp_Pvector(:,j)]=mvSDE(Multi_Data,numSymbol,Depth);
    end
    MorphMatrix_RC=mean(temp_MorphMatrix,3);
    Pvector_RC=mean(temp_Pvector,2);

    apen=0;
    A=[];B=[];
    for k=1:numSymbol^Depth
        for j=1:numSymbol
            if MorphMatrix_RC(k,j)~=0
                A=[A,Pvector_RC(k)*MorphMatrix_RC(k,j)];
                B=[B,MorphMatrix_RC(k,j)];
            end
        end
    end
    apen=apen-sum(A.*log(A));
    % normalizing the computed entropy
    Out_mvMSDE(i)=apen/log(numSymbol^(Depth+1)); 
end
end

function M_Data = Multi(Data,Scale)
%  generate the consecutive coarse-grained time series
%  Input:   data: time series;
%           scale: the scale factor;
%  Output: 
%           M_Data: the coarse-grained time series at the scale factor S.
%  code is arranged by yyt in 2018.08   yangyuantaohit@163.com
L = size(Data,2);
J = fix(L/Scale);
M_Data=NaN*ones(size(Data,1),J);
for j=1:size(Data,1)
    for i=1:J
        M_Data(j,i) = mean(Data(j,(i-1)*Scale+1:i*Scale));
    end
end
end

function [Out_mvSDE,MorphMatrix_CH,Pvector_CH]=mvSDE(Data,numSymbol,Depth)
%  Inputs:   
%           Data: the time series, a matrix;
%           numSymbol:  the number of the symbol kinds
%           Depth: the depth of the state;
%  Output: 
%           Out_mvMSDE: the Refine Composite SpatioTemporal Entropy
%  code is arranged by yyt in 2018.08   yangyuantaohit@163.com
CH=size(Data,1); % Channel of signal
MorphMatrix=NaN*ones(numSymbol^Depth,numSymbol,CH);
Pvector=NaN*ones(numSymbol^Depth,CH);
for i=1:CH 
    [MorphMatrix(:,:,i),Pvector(:,i)]=SymbolicAnalysis(Data(i,:),numSymbol,Depth);
    % normalize the computed vector
    Pvector(:,i) = Pvector(:,i)/sum(Pvector(:,i))/CH;
    % normalize each row of Matrix to make it a stochastic matrix
    for j=1:numSymbol^Depth
        rowSum = sum(MorphMatrix(j,:,i));
        if(rowSum~=0)
            MorphMatrix(j,:,i) = MorphMatrix(j,:,i)/rowSum/CH;     
        end
    end
end
MorphMatrix_CH=sum(MorphMatrix,3);
Pvector_CH=sum(Pvector,2);

apen=0;
A=[];B=[];
for i=1:numSymbol^Depth
    for j=1:numSymbol
       if MorphMatrix_CH(i,j)~=0
           A=[A,Pvector_CH(i)*MorphMatrix_CH(i,j)];
           B=[B,MorphMatrix_CH(i,j)];
       end
    end
end 
apen=apen-sum(A.*log(A));
% normalizing the computed entropy
Out_mvSDE=apen/log(numSymbol^(Depth+1)); 
end



function [MorphMatrix,Pvector]=SymbolicAnalysis(Data,numSymbol,Depth)
% Compute the Probability Vector
% using normal cumulative distribution function (NCDF)
% Input:   
%         Data1:        the raw data - a vector of size 1 x N (the number of sample points)
%         Data2:        the raw data - a vector of size 1 x N (the number of sample points)
%         Depth:       the embedding demension
%         numSymbol:   the number of symbol
% Output:   
%         PVec:        the Probability Vector of data
%         MorphMatrix: the Morph Matrix of data
%  code is arranged by yyt in 2018.08   yangyuantaohit@163.com
if min(size(Data))~=1
    Data=Data(:);
end

L = length(Data);
y = normcdf(Data,mean(Data),std(Data));

epsilon=1e-10;
y(y>1-epsilon)=1-epsilon;
y(y<epsilon)=epsilon;
% for i=1:L
%     if y(i)>1-epsilon
%         y(i)=1-epsilon;
%     end
%     if y(i)<epsilon
%         y(i)=epsilon;
%     end
% end
symbolSeq=round(y*numSymbol+0.5);
%% generate the embedding vector and compute the pVector
numStates=numSymbol^Depth;        
MorphMatrix = zeros(numStates,numSymbol);
Pvector = zeros(1,numStates);
for jj = 1:L-Depth
    WindowSet=symbolSeq(jj:jj+Depth);
    RowNum=1;
    for i=1:Depth
        RowNum=RowNum+(WindowSet(i)-1)*numSymbol^(Depth-i);
    end
    ColumnNum=WindowSet(1+Depth);
    Pvector(RowNum)=1+Pvector(RowNum);
    MorphMatrix(RowNum,ColumnNum)=1+MorphMatrix(RowNum,ColumnNum);   
end
end
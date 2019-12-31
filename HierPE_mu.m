function HPE=HierPE_mu(data,m,t,n)
%  Calculate the Hierarchical Permutation Entropy (HPE)
%  Input:   data: time series;
%           m: embedding dimension;
%           t: time delay of permuation entropy;
%           n: floor of the Hierarchical Entropy;
%  Output: 
%           HPE: hierarchical permutation entropy.
%  code is arranged by yyt in 2016.07     yangyuantaohit@163.com
if size(data,1)==1
   data=data(:);
end
HPE =zeros(1,2^n);
for i=0:(2^n-1)
    datani = HierarchicalEn(n,i,data);
    HPE(i+1) = PermutationEntropy(datani,m,t);
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

function [apen] = PermutationEntropy(data,m,t)
%  Calculate the permutation entropy
%  Input:   data: time series;
%           m: order of permuation entropy;
%           t: delay time of permuation entropy; 
%  Output:  
%           apen: permuation entropy.
%Ref: G Ouyang, J Li, X Liu, X Li, Dynamic Characteristics of Absence EEG Recordings with Multiscale Permutation %  
%                             Entropy Analysis, Epilepsy Research, doi: 10.1016/j.eplepsyres.2012.11.003
%     X Li, G Ouyang, D Richards, Predictability analysis of absence seizures with permutation entropy, Epilepsy %  
%                            Research,  Vol. 77pp. 70-74, 2007

% 列向量转为行向量
[~,a]=size(data);
if a==1
   data=data';
end 

N = length(data);
permlist = perms(1:m);
c(1:length(permlist))=0;
    
 for i=1:N-t*(m-1)
     [~,iv]=sort(data(i:t:i+t*(m-1)));
     for jj=1:length(permlist)
         if (abs(permlist(jj,:)-iv))==0
             c(jj) = c(jj) + var(data(i:t:i+t*(m-1)),1) ;
         end
     end
 end

hist = c; 
c=hist(hist~=0);
p = c/sum(c);
pe = -sum(p .* log(p));
% normalized
apen=pe/log(factorial(m));
end
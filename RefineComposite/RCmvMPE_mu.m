function [Out_mvMPE]=RCmvMPE_mu(X,m,Scale)
%
% To improve the stabily and reliabilty of multivaraite multiscale entropy (mvMSE) based on mean (mvMSE_mu), especially for short signals, we proposed refined composite mvMSE (RCmvMSE_mu). In RCmvMSE_mu, for scale
%factor tau, tau different multivaraite time series, corresponding to different starting points of the coarse-graining process are created and the RCmvMSE_mu value is defined based on
%the averages of the total number of m- and m+1- dimensional matched vector pairs of those shifted sequences.
%
%
% Inputs:
% X: multivariate signal - a matrix of size nvar (the number of channels) x nsamp (the number of sample points for each channel)
% m: scalar embedding value

% Scale: the number of scale factors
%
% Output:
% Out_mvMpPE: a scalar quantity

t=1; % time delay is equal to 1
Out_mvMPE(1) = MVPE(X,m,t);
for ii=2:Scale
    temp_pp=[];
    for iii=1:ii
        Xs = Multi(X(:,iii:end),ii);
        [e,pp] = MVPE(Xs,m,t);
        temp_pp=[temp_pp; pp];
    end
    
f_matrix=temp_pp; 
%     size(temp_ppp)
    AA=mean(temp_pp,1);
    BB=AA(AA~=0);  
    Out_mvMPE(ii)=-sum(BB .* log(BB))/log(factorial(m));
end
end





function M_Data = Multi(Data,S)

%  generate the consecutive coarse-grained time series based on mean
%  Input:   Data: time series;
%           S: the scale factor

% Output:
%           M_Data: the coarse-grained time series at the scale factor S

L = size(Data,2);
J = fix(L/S);

for j=1:size(Data,1)
    for i=1:J
        M_Data(j,i) = mean(Data(j,(i-1)*S+1:i*S));
    end
end
end




function [Out_MVPE,pp]=MVPE(x,m,t)

CH=size(x,1);
ly=size(x,2); % Length of each signal

for i=1:CH   
 y=x(i,:);  
permlist = perms(1:m);
c(1:length(permlist))=0;
    
 for j=1:ly-t*(m-1)
     [a,iv]=sort(y(j:t:j+t*(m-1)));
     for jj=1:length(permlist)
         if (abs(permlist(jj,:)-iv))==0
             c(jj) = c(jj) + var(y(j:t:j+t*(m-1)),1) ;
         end
     end
 end

hist = c;
 
p(i,:) = c/(sum(c)*CH);%每一个除以总的概率
end

pp=sum(p,1);
ppp=pp(pp~=0);  
Out_MVPE=-sum(ppp .* log(ppp));
Out_MVPE=Out_MVPE/log(factorial(m));

end
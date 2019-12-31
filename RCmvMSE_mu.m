function Out_mvMSE=RCmvMSE_mu(X,sm,sr,stau,Scale)
%
% To improve the stabily and reliabilty of multivaraite multiscale entropy (mvMSE) based on mean (mvMSE_mu), especially for short signals, we proposed refined composite mvMSE (RCmvMSE_mu). In RCmvMSE_mu, for scale
%factor tau, tau different multivaraite time series, corresponding to different starting points of the coarse-graining process are created and the RCmvMSE_mu value is defined based on
%the averages of the total number of m- and m+1- dimensional matched vector pairs of those shifted sequences.
%
%
% Inputs:
% X: multivariate signal - a matrix of size nvar (the number of channels) x nsamp (the number of sample points for each channel)
% sm: scalar embedding value
% sr: scalar threshold value (it is usually equal to 0.15)
% stau: scalar time lag  value (it is usually equal to 1)
% Scale: the number of scale factors
%
% Output:
% Out_mvMSE: a scalar quantity
%
% The code is based on the publicly-available code in
% "http://www.commsp.ee.ic.ac.uk/~mandic/research/Complexity_Stuff.htm" prepared by Prof. Mandic's group
%
% Ref:
% [1] H. Azami and J. Escudero, "Refined Composite Multivariate Generalized Multiscale Fuzzy Entropy:
% A Tool for Complexity Analysis of Multichannel Signals", Physica A, 2016.
% [2] M. U. Ahmed and D. P. Mandic, "Multivariate multiscale entropy
% analysis", IEEE Signal Processing Letters, vol. 19, no. 2, pp.91-94.2012
%
% If you use the code, please make sure that you cite references [1] and [2].
%
% Hamed Azami and Javier Escudero Rodriguez
% hamed.azami@ed.ac.uk and javier.escudero@ed.ac.uk
%
%  10-June-16



%%
% Because multi-channel signals may have different amplitude ranges, the distances calculated on embedded vectors may be biased
%toward the largest amplitude ranges variates. For this reason, we scale all the data channels to the same amplitude range and
% we normalize each data channel to unit standard deviation so that the total variation becomes equal to the number of channels
% or variables [1],[2].

X=zscore(X');
r = sr*(sum(std(X)));
M = sm*ones(1,size(X,2)); % M: embedding vector
tau=ones(1,size(X,2))*stau; % tau: time lag vector
X=X';

Out_mvMSE=mvSE(X,M,r,tau);


for ii=2:Scale
    temp_A=[];
    temp_B=[];
    for iii=1:ii
        Xs = Multi(X(:,iii:end),ii);
        [e,A,B] = mvSE(Xs,M,r,tau)   ;
        temp_A=[temp_A A];
        temp_B=[temp_B B];
    end
    AA=sum(temp_A);
    BB=sum(temp_B);
    Out_mvMSE(ii)=log(AA/BB);
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
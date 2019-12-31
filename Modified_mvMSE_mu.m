function Out_mvMSE=Modified_mvMSE_mu(X,sm,sr,stau,Scale)
% This function calculates multivariate multiscale (sample) entropy (mvMSE) whose
% coarse-graining is based on mean
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

Out_mvMSE(1)=mvSE(X,M,r,tau);
for k=2:Scale
    Xs = Multi(X,k);
    Out_mvMSE(k) = mvSE(Xs,M,r,tau);
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
        M_Data(j,i) = mean(Data(j,(i:i+S-1)));
    end
end

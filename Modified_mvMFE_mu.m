function Out_mvMFE=Modified_mvMFE_mu(X,sm,sr,n,stau,Scale)
% This function calculates multivariate multiscale fuzzy entropy (mvMFE) whose
% coarse-graining is based on mean
%
%
% Inputs:
% X: multivariate signal - a matrix of size nvar (the number of channels) x nsamp (the number of sample points for each channel)
% sm: scalar embedding value
% sr: scalar threshold value (it is usually equal to 0.15)
% n: fuzzy power (it is usually equal to 2)
% stau: scalar time lag  value (it is usually equal to 1)
% Scale: the number of scale factors
%
% Output:
% Out_mvMFE: a scalar quantity
%
% Ref:
% [1] H. Azami and J. Escudero, "Refined Composite Multivariate Generalized Multiscale Fuzzy Entropy:
% A Tool for Complexity Analysis of Multichannel Signals", Physica A, 2016.
%
% If you use the code, please make sure that you cite reference [1].
%
% Hamed Azami and Javier Escudero Rodriguez
% hamed.azami@ed.ac.uk and javier.escudero@ed.ac.uk
%
%  10-June-16



%%
% Because multi-channel signals may have different amplitude ranges, the distances calculated on embedded vectors may be biased
%toward the largest amplitude ranges variates. For this reason, we scale all the data channels to the same amplitude range and
% we normalize each data channel to unit standard deviation so that the total variation becomes equal to the number of channels
% or variables [1].

X=zscore(X');
r = sr*(sum(std(X)));
M = sm*ones(1,size(X,2)); % M: embedding vector
tau=ones(1,size(X,2))*stau; % tau: time lag vector
X=X';

Out_mvMFE(1)=mvFE(X,M,r,n,tau);
for k=2:Scale
    Xs = Multi(X,k);
    Out_mvMFE(k) = mvFE(Xs,M,r,n,tau);
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
        M_Data(j,i) = mean(Data(j,(i:i+S-1)));
    end
end
end
function Out_RCMDE=RCMDE(x,m,c,tau,Scale)
%
% This function calculates the refined composite multiscale dispersion entropy (RCMDE) of a univariate signal x
%
% Inputs:
%
% x: univariate signal - a vector of size 1 x N (the number of sample points)
% m: embedding dimension
% c: number of classes (it is usually equal to a number between 3 and 9 - we used c=6 in our studies)
% tau: time lag (it is usually equal to 1)
% Scale: number of scale factors
%
%Outputs:
%
% Out_RCMDE: a vector of size 1 * Scale - the RCMDE of x
%
%
% Ref:
% [1] H. Azami, M. Rostaghi, D. Abasolo, and J. Escudero, "Refined Composite Multiscale Dispersion Entropy and its Application to Biomedical
% Signals", IEEE Transactions on Biomedical Engineering, 2017.
% [2] M. Rostaghi and H. Azami, "Dispersion Entropy: A Measure for Time-Series Analysis", IEEE Signal Processing Letters. vol. 23, n. 5, pp. 610-614, 2016.
%
% If you use the code, please make sure that you cite references [1] and [2].
%
% Hamed Azami and Javier Escudero Rodriguez
% Emails: hamed.azami@ed.ac.uk and javier.escudero@ed.ac.uk
%
%  20-January-2017
%%

Out_RCMDE=NaN*ones(1,Scale);

Out_RCMDE(1)=DisEn_NCDF(x,m,c,tau);

sigma=std(x);
mu=mean(x);

for j=2:Scale
    pdf=[];
    for jj=1:j
        xs = Multi(x(jj:end),j);
        [DE, T_pdf]=DisEn_NCDF_ms(xs,m,c,mu,sigma,tau);
        pdf=[pdf ; T_pdf];
    end
    pdf=mean(pdf,1);
    pdf=pdf(pdf~=0);
    Out_RCMDE(j)=-sum(pdf .* log(pdf));
end
end

function M_Data = Multi(Data,S)

%  generate the consecutive coarse-grained time series
%  Input:   Data: time series;
%           S: the scale factor
% Output:
%           M_Data: the coarse-grained time series at the scale factor S

L = length(Data);
J = fix(L/S);

for i=1:J
    M_Data(i) = mean(Data((i-1)*S+1:i*S));
end
end

function [Out_DisEn, npdf]=DisEn_NCDF(x,m,nc,tau)
%
% This function calculates dispersion entropy (DisEn) of a univariate
% signal x, using normal cumulative distribution function (NCDF)
%
% Inputs:
%
% x: univariate signal - a vector of size 1 x N (the number of sample points)
% m: embedding dimension
% nc: number of classes (it is usually equal to a number between 3 and 9 - we used c=6 in our studies)
% tau: time lag (it is usually equal to 1)
%
% Outputs:
%
% Out_DisEn: scalar quantity - the DisEn of x
% npdf: a vector of length nc^m, showing the normalized number of disersion patterns of x
%
% Ref:
%
% [1] H. Azami, M. Rostaghi, D. Abasolo, and J. Escudero, "Refined Composite Multiscale Dispersion Entropy and its Application to Biomedical
% Signals", IEEE Transactions on Biomedical Engineering, 2017.
% [2] M. Rostaghi and H. Azami, "Dispersion Entropy: A Measure for Time-Series Analysis", IEEE Signal Processing Letters. vol. 23, n. 5, pp. 610-614, 2016.
%
% If you use the code, please make sure that you cite references [1] and [2].
%
% Hamed Azami, Mostafa Rostaghi, and Javier Escudero Rodriguez
% hamed.azami@ed.ac.uk, rostaghi@yahoo.com, and javier.escudero@ed.ac.uk
%
%  20-January-2017
%%


N=length(x);

sigma_x=std(x);
mu_x=mean(x);

y=normcdf(x,mu_x,sigma_x);

for i_N=1:N
    if y(i_N)==1
        y(i_N)=1-(1e-10);
    end
    
     if y(i_N)==0
        y(i_N)=(1e-10);
     end
end

z=round(y*nc+0.5);

all_patterns=[1:nc]';

for f=2:m
    temp=all_patterns;
    all_patterns=[];
    j=1;
    for w=1:nc
        [a,b]=size(temp);
        all_patterns(j:j+a-1,:)=[temp,w*ones(a,1)];
        j=j+a;
    end
end

for i=1:nc^m
    key(i)=0;
    for ii=1:m
        key(i)=key(i)*10+all_patterns(i,ii);
    end
end


embd2=zeros(N-(m-1)*tau,1);
for i = 1:m, 
    embd2=[z(1+(i-1)*tau:N-(m-i)*tau)]'*10^(m-i)+embd2;
end


pdf=zeros(1,nc^m);

for id=1:nc^m
    [R,C]=find(embd2==key(id));
    pdf(id)=length(R);
end

npdf=pdf/(N-(m-1)*tau);
p=npdf(npdf~=0);
Out_DisEn = -sum(p .* log(p));
end

function [Out_DisEn, npdf]=DisEn_NCDF_ms(x,m,nc,mu,sigma,tau)
%
% This function calculates dispersion entropy (DisEn) using normal cumulative distribution function (NCDF) with defined mean (mu) and standard deviation (sigma) values.
%
% Inputs:
%
% x: univariate signal - a vector of size 1 x N (the number of sample points)
% m: embedding dimension
% nc: number of classes (it is usually equal to a number between 3 and 9 - we used c=6 in our studies)
% tau: time lag (it is usually equal to 1)
%
% Outputs:
%
% Out_DisEn: scalar quantity - the DisEn of x
% npdf: a vector of length nc^m, showing the normalized number of disersion patterns of x
%
% Ref:
%
% [1] H. Azami, M. Rostaghi, D. Abasolo, and J. Escudero, "Refined Composite Multiscale Dispersion Entropy and its Application to Biomedical
% Signals", IEEE Transactions on Biomedical Engineering, 2017.
% [2] M. Rostaghi and H. Azami, "Dispersion Entropy: A Measure for Time-Series Analysis", IEEE Signal Processing Letters. vol. 23, n. 5, pp. 610-614, 2016.
%
% If you use the code, please make sure that you cite references [1] and [2].
%
% Hamed Azami, Mostafa Rostaghi, and Javier Escudero Rodriguez
% hamed.azami@ed.ac.uk, rostaghi@yahoo.com, and javier.escudero@ed.ac.uk
%
%  20-January-2017
%%
N=length(x);

y=normcdf(x,mu,sigma);

for i_N=1:N
    if y(i_N)==1
        y(i_N)=1-(1e-10);
    end
    
     if y(i_N)==0
        y(i_N)=(1e-10);
     end
end

z=round(y*nc+0.5);

all_patterns=[1:nc]';

for f=2:m
    temp=all_patterns;
    all_patterns=[];
    j=1;
    for w=1:nc
        [a,b]=size(temp);
        all_patterns(j:j+a-1,:)=[temp,w*ones(a,1)];
        j=j+a;
    end
end

for i=1:nc^m
    key(i)=0;
    for ii=1:m
        key(i)=key(i)*10+all_patterns(i,ii);
    end
end

embd2=zeros(N-(m-1)*tau,1);
for i = 1:m, 
    embd2=[z(1+(i-1)*tau:N-(m-i)*tau)]'*10^(m-i)+embd2;
end

pdf=zeros(1,nc^m);

for id=1:nc^m
    [R,C]=find(embd2==key(id));
    pdf(id)=length(R);
end

npdf=pdf/(N-(m-1)*tau);
p=npdf(npdf~=0);
Out_DisEn = -sum(p .* log(p));
end


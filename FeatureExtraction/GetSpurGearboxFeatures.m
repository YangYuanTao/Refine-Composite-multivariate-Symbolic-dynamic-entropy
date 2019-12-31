function Feature_Paras=GetSpurGearboxFeatures(Vibdata,Fs,fn,ratio,Z)
% Get_Features4Gearbox       Function: Compute the features for gearbox
% fault diagnosis
% Inpute：
%      （1）Vibdata —— Vibration signal  organized in column;
%      （2）Fs —— Sampling frequency (Hz);
%      （3）fn —— the rotation frequency of input shaft;
%      （4）ratio —— Gear reduction ratio from input shaft to idler shaft; 
%      （5）Z —— the number of teeth ;

% Output：
%      （1）Feature_Paras—— the corresponding feature parameters to Vibdata
%      （2）RelSig—— related signals corresponding to Vibdata

% References:
% [ Samuel P D, Pines D J. A review of vibration-based techniques for 
%   helicopter transmission diagnostics[J]. Journal of sound and vibration,
%   2005, 282(1): 475-508.].

% [Cerrada M, Zurita G, Cabrera D, et al. Fault diagnosis in spur gears 
%  based on genetic algorithm and random forest[J]. Mechanical Systems and 
%  Signal Processing, 2016, 70: 87-103.]

% [	Tom, K.F., Survey of Diagnostic Techniques for Dynamic Components.
%       Technical Report of Army Research Laboratory, 2010.]

% Zheng Huailiang 2017-09-10
% Yang Yuantao 20191221
% v2 ------ the features are computed using original vibaration signal
% v4 ------ 在计算加速度信号的特征前，去掉其均值； 2018/12/27

if nargin < 5
    error('Not enough input!');
end

% Note: in the following discuss, only one harmonic of fn and fs are
% considered

% compute meshing frequency
fn_i = fn*ratio;        % rotation frequency of idler shaft
fs = fn_i*Z;            % meshing frequency    

L=length(Vibdata);  % the number of samples
Vibdata = bsxfun(@minus, Vibdata, mean(Vibdata,1));

%% ===== Special features for gear fault diagnosis  =====%%
%% ========1. 13 time domain feature======= %%
%计算有量纲统计特征值：峰值（Peak）、峰峰值（PPV）、均值（MV）、
%平均幅值（MA）、方根幅值（SRA）、标准差（SD）、有效值（均方根值）（RMS）
Peak=max(abs(Vibdata));
PPV=max(Vibdata)-min(Vibdata);
MV=mean(Vibdata);
MA=mean(abs(Vibdata));
SRA=mean(sqrt(abs(Vibdata)))^2;
SD=std(Vibdata,1);
RMS=sqrt(mean(Vibdata.^2));
%计算无量纲统计特征值：偏度指标（SKE）、峭度指标(KUR)、峰值指标(CRE)、脉冲指标(IMP)、裕度指标(CLE)、波形指标（BY）
SKE=mean((Vibdata-MV).^3)/RMS^3;
KUR=mean((Vibdata-MV).^4)/RMS^4;
CRE=Peak/RMS;
IMP=Peak/MA;
CLE=Peak/SRA;
BY=RMS/MA;
% construct time domain feature
time_paras=NaN*ones(13,1);
time_paras(1,1)=Peak; time_paras(2,1)=PPV;time_paras(3,1)=MV;
time_paras(4,1)=MA; time_paras(5,1)=SRA;time_paras(6,1)=SD;
time_paras(7,1)=RMS; time_paras(8,1)=SKE;time_paras(9,1)=KUR;
time_paras(10,1)=CRE; time_paras(11,1)=IMP;time_paras(12,1)=CLE; 
time_paras(13,1)=BY; 

% compute the regular component signal and difference signal
% !!!! attention the sample rate for tsa and roration frquency of idler shaft !!!
[d1,d] = RelatedSignal_d(Vibdata,Fs,fs,fn_i); 
% compute the re signal
re = RelatedSignal_re(Vibdata);

%==========2. FM4============%
mean_d = mean(d);
minus_d = bsxfun(@minus, d, mean_d);
FM4 = L*sum(minus_d.^4)./(sum(minus_d.^2).^2);

%==========3. M6A============%
M6A = L*L*sum(minus_d.^6)./(sum(minus_d.^2).^3);

%=========4. Energy ratio, ER========%
ER = sqrt(mean(d.^2)./mean(d1.^2));

%=========5. Energy Operator, EOP========%
mean_re = mean(re);
minus_re = bsxfun(@minus, re, mean_re);
EOP = L*sum(minus_re.^4)./(sum(minus_re.^2).^2);

%% ===== 6. Features in frequency domain =====%%
% the frequency domain features are computed based on orignal vibarationsignal
f=Fs/2*linspace(0,1,round(L/2)+1);  % the frequency axis(row vector)
YF = fft(Vibdata,L)/L;
FZYF=2*abs(YF(1:round(L/2)+1));     % fft
% the statistical features in frequency domain
FC = (f*FZYF)/sum(FZYF);                 % Frequency centre
RMSF = sqrt((f.^2*FZYF)/sum(FZYF));      % Root mean square frequency
STDF = sqrt(((f-FC).^2*FZYF)/sum(FZYF)); % Standard deviation frequency
CP1 = mean(FZYF'.*(f-FC).^3)/(STDF^3);
CP2 = STDF/FC;
CP3 = mean(FZYF'.*(f-FC).^4)/(STDF^4);

fre_paras = NaN*ones(10,1);
fre_paras(1,1) = FC; fre_paras(2,1) = RMSF; fre_paras(3,1) = STDF;
fre_paras(4,1) = CP1; fre_paras(5,1) = CP2; fre_paras(6,1) = CP3;

% the 1st rotation frequency amplitude of input shaft
Am_fn = NaN*ones(3,1);
eps1=6;                              %设置一阶转频容许变化精度
n1 = fix((fn-eps1)*L/Fs)+1;          %得到找实际一阶转频的容许区间
n2 = ceil((fn+eps1)*L/Fs)+1;
[Am_fn(1,1),I] = max(FZYF(n1:n2));   %找到容许区间的最大值并作为实际的一阶转频
IFN = n1+I-1;
% the 2nd 3rd rotation frequency amplitude of input shaft
eps2 = 6;                               %设置多阶转频容许变化精度
for k = 2:3
    n1=fix((k*f(IFN)-eps2)*L/Fs)+1;     %得到找实际多阶转频的容许区间
    n2=ceil((k*f(IFN)+eps2)*L/Fs)+1;
    [Am_fn(k,1),~]=max(FZYF(n1:n2));    %找到容许区间的最大值并作为实际的多阶转频
end
fre_paras(7:9,1) = Am_fn; 
% the meshing frequency amplitude
eps1=8; 
n1 = fix((fs-eps1)*L/Fs)+1;
n2 = ceil((fs+eps1)*L/Fs)+1;
[fre_paras(10,1),~] = max(FZYF(n1:n2));

% store the feature prarmeter
Feature_Paras=NaN*ones(10,1);
Feature_Paras(1:13,1) = time_paras;
Feature_Paras(14:23,1) = fre_paras;
Feature_Paras(24,1) = FM4;
Feature_Paras(25,1) = M6A;
Feature_Paras(26,1) = ER;
Feature_Paras(27,1) = EOP;

% store the related signal
% RelSig = cell(3,1);
% RelSig{1,1} = d1;
% RelSig{2,1} = d;
% RelSig{3,1} = re;
end


function [d1,d] = RelatedSignal_d(Vibdata,Fs,fs,fn)
% RelatedSignal_d.m       Function: Compute the related signal for
% computing the special parameters of gear damage detection
% Inpute：
%      （1）Vibdata —— Vibration signal  organized in column;
%      （2）Fs —— Sampling frequency (Hz);
%      （3）fs —— the meashing frequency;
%      （4）fn —— the rotation frequency;
% Output：
%      （1）s —— the envelope signal of band-passed filtered signal about
%      meashing frequency
%      （2）r —— Residual signal
%      （3）d1 —— Regular meshing components signal
%      （4）d —— Difference signal
%      （5）re —— x(i)^2-x(i-1)*x(i+1) signal

if nargin < 4
    error('The wrong input!')
end

[N,Num] = size(Vibdata);

eps1 = 3; % to 3.0H/F  Rp = 3;Rs=5

%%======3. Regular meshing components signal & difference signal ======%%
d1 = zeros(N,Num);
%  construct band-passed filter - store rotation frequency component
wp=[fn-eps1,fn+eps1];
ws=[fn-2*eps1,fn+2*eps1];
Rp=3;
Rs=5;
[n,Wn]=buttord(wp*2/Fs,ws*2/Fs,Rp,Rs,'s');
[b,a]=butter(n,Wn,'bandpass');

for i = 1:Num
    d1(:,i) = filter(b,a,Vibdata(:,i));
end

%  construct band-passed filter - store meshing frequency component and
%  fisrt-order bandside
wp=[fs-fn,fs+fn];
ws=[fs-2*fn,fs+2*fn];
Rp=3;
Rs=10;
[n,Wn]=buttord(wp*2/Fs,ws*2/Fs,Rp,Rs,'s');
[b,a]=butter(n,Wn,'bandpass');

for i = 1:Num
    d1(:,i) = d1(:,i)+filter(b,a,Vibdata(:,i));
end

d = Vibdata - d1;

end


function re = RelatedSignal_re(Vibdata)
% RelatedSignals.m       Function: Compute the related signal for
% computing the special parameters of gear damage detection
% Inpute：
%      （1）Vibdata —— Vibration signal  organized in column;
%      （2）Fs —— Sampling frequency (Hz);
%      （3）fs —— the meashing frequency;
%      （4）fn —— the rotation frequency;
% Output：
%      （1）s —— the envelope signal of band-passed filtered signal about
%      meashing frequency
%      （2）r —— Residual signal
%      （3）d1 —— Regular meshing components signal
%      （4）d —— Difference signal
%      （5）re —— x(i)^2-x(i-1)*x(i+1) signal

if nargin < 1
    error('The wrong input!')
end

N= size(Vibdata,1);

%%======4. Re signal ======%%
% re = zeros(N,Num);
x11 = [Vibdata(end,:);Vibdata;Vibdata(1,:)];
re = x11(2:N+1,:).^2 - x11(1:N,:).*x11(3:N+2,:);
end




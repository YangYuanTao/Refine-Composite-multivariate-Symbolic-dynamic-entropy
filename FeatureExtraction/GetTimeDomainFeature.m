function SYPara=GetTimeDomainFeature(x)
%    GetTimeDomainFeature.m      函数：计算给定信号的时域统计特征值
%           x                    输入:要计算时域特征参数的信号(没有时间列)
%        SYPara                  输出：计算得到的时域统计特征参数（13×1矩阵）
%     SYPara(1,1)-峰值（Peak）; SYPara(2,1)-峰峰值（PPV）;SYPara(3,1)-均值（MV）
%     SYPara(4,1)-平均幅值（MA）; SYPara(5,1)-方根幅值（SRA）;SYPara(6,1)-标准差（SD）
%     SYPara(7,1)-有效值（RMS）; SYPara(8,1)-偏度指标（SKE）;SYPara(9,1)-峭度指标(KUR)
%     SYPara(10,1)-峰值指标(CRE); SYPara(11,1)-脉冲指标(IMP);SYPara(12,1)-裕度指标(CLE)
%     SYPara(13,1)-波形指标(BY)

if nargin~=1
    error('The Wrong Input');
end

%计算有量纲统计特征值：峰值（Peak）、峰峰值（PPV）、均值（MV）、
%平均幅值（MA）、方根幅值（SRA）、标准差（SD）、有效值（均方根值）（RMS）
Peak=max(abs(x));
PPV=max(x)-min(x);
MV=mean(x);
MA=mean(abs(x));
SRA=mean(sqrt(abs(x)))^2;
SD=std(x,1);
RMS=sqrt(mean(x.^2));

%计算无量纲统计特征值：偏度指标（SKE）、峭度指标(KUR)、峰值指标(CRE)、脉冲指标(IMP)、裕度指标(CLE)、波形指标（BY）
SKE=mean((x-MV).^3)/RMS^3;
KUR=mean((x-MV).^4)/RMS^4;
CRE=Peak/RMS;
IMP=Peak/MA;
CLE=Peak/SRA;
BY=RMS/MA;
    
SYPara=zeros(12,1);
SYPara(1,1)=Peak; SYPara(2,1)=PPV;SYPara(3,1)=MV;
SYPara(4,1)=MA; SYPara(5,1)=SRA;SYPara(6,1)=SD;
SYPara(7,1)=RMS; SYPara(8,1)=SKE;SYPara(9,1)=KUR;
SYPara(10,1)=CRE; SYPara(11,1)=IMP;SYPara(12,1)=CLE; 
SYPara(13,1)=BY; 

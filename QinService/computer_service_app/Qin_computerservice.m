function [StructDate,UI_printInf,flag_Qin,flag_isUsed]=Qin_computerservice(Sigframe_decoded,Pf,flag_yundong)
%  flag_isUsed 结构体数据中存储在StructDate{1,11}的数据结构体是否有效
%  StructDate  结构体数据  元胞矩阵 1*11

% 添加路径依赖项
addpath 'E:\A_Matlab2020a\Matlab2020a\bin\computer_service_app'


% 场景选择 优先读取Q/T数据
flag_Qin1=0;       % 情景标志位
flag_Qin2=0;
flag_Qin3=0; 
count_zhen = 1;    % 暂时默认只处理一帧
Q_T=zeros(2,count_zhen);  
sitution=zeros(1,count_zhen);      % 用于说明当前帧适应的场景类型
for i=1:count_zhen
    
    Q_T(1,i)=Sigframe_decoded{1,20};   % 获取Q散射体数量
    Q_T(2,i)=Sigframe_decoded{1,56};   % 获取T散射体数量
    if Q_T(1,i)==0 && Q_T(2,i)==0
        flag_Qin1=flag_Qin1+1;
        sitution(i)=1;
    elseif Q_T(1,i)> 0 && Q_T(2,i)==0
         flag_Qin2=flag_Qin2+1;
        sitution(i)=2;
    elseif Q_T(1,i)> 0 && Q_T(2,i)> 0
          flag_Qin3=flag_Qin3+1;
        sitution(i)=3;
    end
end
if flag_Qin1>(count_zhen/2)                % 确定定位场景 针对多帧Q-T值不统一时选择定位场景最多
    flag_Qin1=1; flag_Qin2=0;flag_Qin3=0;
    flag_Qin =1 ;                          % UI获知定位情景
elseif flag_Qin2>(count_zhen/2)
    flag_Qin1=0; flag_Qin2=1;flag_Qin3=0;
    flag_Qin =2;
elseif flag_Qin3>(count_zhen/2)
    flag_Qin1=0; flag_Qin2=0;flag_Qin3=1;
    flag_Qin =3;
else
    flag_Qin1=0; flag_Qin2=1;flag_Qin3=0;  %% 无优势帧选择情景2
    flag_Qin =2;
end
% 进入情景计算
% 函数返回参数说明： 
if flag_Qin1==1
    disp('情景1状态Mention：进入到升空散射体定位场景1中');
    
elseif flag_Qin2==1
     % 情景2计算
    disp('情景2状态Mention：进入到升空散射体定位场景2中');
    [StructDate,UI_printInf,flag_isUsed]=Qin2_computerservice(Sigframe_decoded,Pf,flag_yundong); %% 还未考虑多帧的情况
    
elseif flag_Qin3==1
     % 情景3计算
     disp('情景3状态Mention：进入到升空散射体定位场景3中');
     [StructDate,UI_printInf,flag_isUsed]=Qin3_computerservice(Sigframe_decoded,Pf,flag_yundong);
end



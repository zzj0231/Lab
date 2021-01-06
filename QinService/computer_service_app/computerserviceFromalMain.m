function [StructDate,UI_printInf,flag_Qin,flag_isUsed] = computerserviceFromalMain(Sigframe_decoded,parameterSit1_cell,parameterSit2_cell,isMoveTarget)
% 计算服务器三种情景应用集合函数

%  parameterSit1_cell： 移动单站纯DOA定位方式
%  parameterSit2_cell： 升空散射体辅助定位方式
%  flag_isUsed 结构体数据中存储在StructDate{1,11}的数据结构体是否有效
%  StructDate  结构体数据  元胞矩阵 1*11
%  UI_printInf 待打印的定位过程信息

% 加路径依赖项
addpath 'E:\A_Matlab2020a\Matlab2020a\bin\computer_service_app'

% 场景选择 优先读取Q/T数据
flag_Qin2=0;
flag_Qin3=0; 
count_zhen = 1;    % 暂时默认只处理一帧
Q_T=zeros(2,count_zhen);  
sitution=zeros(1,count_zhen);      % 用于说明当前帧适应的场景类型

if isMoveTarget==1
    flag_Qin =1;
    disp('情景1状态Mention：进入到升空散射体定位场景1中');
    Num_m = parameterSit1_cell{1};
    Num_zhen = parameterSit1_cell{2};
    ya = parameterSit1_cell{3};
    yb = parameterSit1_cell{4};
    r=parameterSit1_cell{5};
    k_up =parameterSit1_cell{6};
    k_down =parameterSit1_cell{7};
    [StructDate,flag_isUsed]=Qin1_computerseverceFormal(MutilState_cell,Num_m,Num_zhen,ya, yb, k_up, k_down, r);    % 定位结果 当前CPI瞬时结构体 二进制结构图(未定义起始和截止 缺少定位结果字节) 
else
    for i=1:count_zhen
        Q_T(1,i)=Sigframe_decoded{1,20};   % 获取Q散射体数量
        Q_T(2,i)=Sigframe_decoded{1,56};   % 获取T散射体数量
        if Q_T(1,i)> 0 && Q_T(2,i)==0
            flag_Qin2=flag_Qin2+1;
            sitution(i)=2;
        elseif Q_T(1,i)> 0 && Q_T(2,i)> 0
            flag_Qin3=flag_Qin3+1;
            sitution(i)=3;
        end
    end
    % UI获知定位情景
    if flag_Qin2>(count_zhen/2)
        flag_Qin2=1;flag_Qin3=0;
        flag_Qin =2;
    elseif flag_Qin3>(count_zhen/2)
        flag_Qin2=0;flag_Qin3=1;
        flag_Qin =3;
    else
        flag_Qin2=1;flag_Qin3=0;  % 无优势帧选择情景2
        flag_Qin =2;
    end
    
    %进入情景计算
    if flag_Qin2==1
        % 情景2计算
        disp('情景2状态Mention：进入到升空散射体定位场景2中');
        Pf = parameterSit2_cell{1};
        flag_yundong = parameterSit2_cell{2};
        [StructDate,UI_printInf,flag_isUsed]=Qin2_computerservice(Sigframe_decoded,Pf,flag_yundong);              % 还未考虑多帧的情况
        
    elseif flag_Qin3==1
        % 情景3计算
        disp('情景3状态Mention：进入到升空散射体定位场景3中');
        Pf = parameterSit2_cell{1};
        flag_yundong = parameterSit2_cell{2};
        [StructDate,UI_printInf,flag_isUsed]=Qin3_computerservice(Sigframe_decoded,Pf,flag_yundong);
    end
end
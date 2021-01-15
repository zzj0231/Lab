function [StructDate,UI_printInf,flag_isUsed]=Qin2_computerservice(Sigframe_decoded,Pf,flag_yundong)
%%% 情景2 Q>0 T=0
%%% Sigframe_decoder :计时器收集到的有效帧数据解析后的结果 默认是一帧元胞数据
%%% flag_yundong: 低采样率针对的是运动目标 0静止 1：开启 
%%% 输出参数： 
%%%  Loc_dingwei:    输出定位结果
%%%  StructDate:     定位结构
%%%  flag_isUsed=0   定位结果是否有效

%%% 特征滤波 期望信号再确认
  % there? or other?
addpath 'E:\A_Matlab2020a\Matlab2020a\bin\computer_service_app\ComSer2_situtation2'    

% 截胡准备区
% simulation_sigframe = load('E:\混合数据帧2解析\Mixframe2_2020_11_12_10_53_58_650_decode.mat');
% simulation_sigframe =simulation_sigframe.Struct_data;
% simulation_RD = load('E:\混合数据帧2解析\RD_cell_simulation.mat');
% simulation_RD = simulation_RD.RD_cell;
%

% 获取直达波方向 读取升空散射体方位扇区RD数据块
flag_yundong=0;

Scample_rate_1 = 32500000;        % scamping rate_1: 32.5MHZ
Scample_rate_2 = 16250000;        % scamping rate_2: 16.25MHZ
Scample_rate_3 = 8125000;         % scamping rate_3: 8.125MHZ
Scample_rate_4 = 4062500;         % scamping rate_4: 4.0625MHZ
Scample_rate_5 = 203125;          % scamping rate_5: 203125KHZ
StructDate=cell(1,11);
 StructDate{1}=[Sigframe_decoded{5},Sigframe_decoded{6},Sigframe_decoded{7},Sigframe_decoded{8},Sigframe_decoded{9},Sigframe_decoded{10},Sigframe_decoded{11}]; % 年-月-日-时-分-秒-高精度计数
 Order_num =Sigframe_decoded{4};
 if Order_num==0
     scampling=0.203125;  %%信号采样率2013.125ksps
 elseif Order_num==1
     Band=3;            %% 带宽3MHz
     scampling=4.0625;  %% 信号采样率4.0625Msps
 elseif Order_num==2
     Band=5;            %% 带宽5MHz
     scampling=8.125;   %% 信号采样率8.125Msps
 elseif Order_num==3
     Band=10;           %% 带宽10MHz
     scampling=16.25;   %% 信号采样率16.25Msps
 elseif Order_num==4
     Band=30;           %% 带宽30MHz
     scampling=32.5;    %% 信号采样率32.25Msps
 end
 StructDate{2}=[Sigframe_decoded{3},Band,scampling];  % 中心频率-带宽-采样率
 if scampling ==32.5
     L = 163; K=1;
 elseif scampling ==16.25
     L = 82; K=2;
 elseif scampling ==8.125
     L = 41; K=3;
 elseif scampling ==4.0625
     L = 21; K=6;
 elseif scampling ==0.203125
     L = 1;  K=113;
 end
 StructDate{3}=[1,L,scampling]; % 距离单元序号范围-距离分辨率(ns)
 StructDate{4}=[1,K,6];         % 多普勒单元序号范围-多普勒分辨频率(HZ)
 StructDate{5}=[0,1];           % 工作场景分类
 StructDate{6}=zeros(2,6);      % UCA经度 纬度 高度 速度x y z
                                % 磁x 磁y 磁z 俯仰 横向 0
 UCA_state =[Sigframe_decoded{12},Sigframe_decoded{13},Sigframe_decoded{14}];
 UCA_magnetic = [Sigframe_decoded{17},Sigframe_decoded{18},Sigframe_decoded{19}];
 [UCA_azmiuth,UCA_pitch] = AntennaCor2StandardStationCor(UCA_magnetic,Sigframe_decoded{15},Sigframe_decoded{16},0,0);  % UCA方位角和俯仰角转为站心坐标系下角度值
 StructDate{6}(1,1:3)=UCA_state;
 StructDate{6}(1,1:3)=UCA_magnetic;
 StructDate{6}(1,4:5)=[UCA_azmiuth,UCA_pitch];
 Q_num = Sigframe_decoded{20};
 StructDate{7} = zeros(5,8);      % S经度 S纬度 S高度  S_速度x S_速度y S_速度z 方向角 俯仰角
 StructDate{7}(1,1) = Q_num;
 index =1;
 for i=1:Q_num
     StructDate{7}(i+1,1:3) = [Sigframe_decoded{20+index},Sigframe_decoded{20+index+1},Sigframe_decoded{20+index+2}];
     StructDate{7}(i+1,4) = Sigframe_decoded{20+index+3};            % 方向
     StructDate{7}(i+1,7) = Sigframe_decoded{20+index+4};            % 速度
     index =index+5;
 end
P_num = Sigframe_decoded{41};
StructDate{8} = zeros(8,3);     % DOA信息 8*3 第一行 DOA的数量
                                % 方向角 俯仰角 来波方向强度
StructDate{8}(1,1)=P_num;
index = 1;
for i=1:P_num
    StructDate{8}(i+1,1:3) = [Sigframe_decoded{41+index},Sigframe_decoded{41+index+1},0];
    index =index+2;
end
StructDate{9}=zeros(5,5);                   % zeros 5*5  T散射体信息
refdoa_order = Sigframe_decoded{65};
StructDate{10}=zeros(1,3);                  % 参考通道信号来波方向 DOA信息中的字段
StructDate{10}(1:2) = StructDate{8}(refdoa_order+1,1:2);
StructDate{10}(3) =  refdoa_order;
TDOA_FDOA_Az_DOA=StructDate;                % TDOA-DOA-Az-FDOA数据结构体 AZ参考通道的序号
doa_ref=StructDate{1,10}(1:2);              % 参考通道DOA
%

%%%% 存储TDOA-DOA-Azimuth-DOA结构体
    for i=10:-1:5
        TDOA_FDOA_Az_DOA{i+1}=TDOA_FDOA_Az_DOA{i};             %%空出胞单元5 填充方位扇区分布信息
    end
        TDOA_FDOA_Az_DOA{5}=[1,12,30];                         %% 方位扇区的序号范围、方位分辨率
%%%

%%% 获取IQ信号采样率和带宽
    Caiyang = StructDate{1,2}(1,3);
    if Caiyang<=33.5 && Caiyang >=4.0525
        flag_caiyang=1;
    else
        flag_caiyang=0;
    end
    
    if Caiyang==32.5
        flag_caiyangl=1;
    elseif Caiyang==16.25
        flag_caiyangl=2;
    elseif Caiyang==8.125
        flag_caiyangl=3;
    elseif Caiyang==4.0625
        flag_caiyangl=4;
    elseif Caiyang==0.203125
        flag_caiyangl=0;
    end
%%%
   %%% 获取散射体坐标、获取散射体速度
    Q_array=StructDate{1,7}; % 提取出散射体的数量 经/纬/高 瞬时速度
    cen_state=StructDate{6}(1,1:3);
   % true area    
    S_xyz=LBM_sXYZ((Q_array(2:Q_num+1,1:3))',cen_state);
    S_xyz = Q_array(2:Q_num+1,1:3)';               %  自制的解析帧散射体数据是站心坐标系 不需要坐标转换 
   % temp area %%%
     % S_xyz=StructDate{1,7}(2:Num_Q+1,1:3)';
   %
    vsi_array=Q_array(2:Q_num+1,4);
    fai_array=Q_array(2:Q_num+1,7);
    UCA_sudu =[0 0];                                                 % UCA速度未知     
    S_vxyz=Count_Sv(S_xyz,vsi_array,fai_array,UCA_sudu);             % 计算散射体数据 散射体坐标、散射体瞬时速度、航向角、UCA速度1*2(暂无) 
    StructDate{1,7}(2:Q_num+1,4:6)=S_vxyz';                          % 散射体速度赋给结构体
    TDOA_FDOA_Az_DOA{1,8}(2:Q_num+1,4:6)=S_vxyz';
    RD_data =cell(1,12);                                             % Sigframe_decoded 提取RD谱数据
    for i=68:79
        RD_data{i-67} = Sigframe_decoded{i};
    end
   
   UI_printInf =cell(1,5);
%  第二部分 获取散射体数量 散射体对应方位角/俯仰角 填装Rd_cell CFAR检测 峰值检测   
   Rd_cell=cell(1,4);       % 声明RD数据容器  用于装载每个散射体对应的方位扇区的RD数据
   Pefeng_cell=cell(1,4);   % Pefeng_cell容器 用于装载有效的谱峰搜索结果
   Q=Q_num;                 % 用于记录RD数据谱有效的散射体数量
                                
   record_cfarout = cell(1,Q_num);
   record_noisetype =cell(1,Q_num);

   count_avail=0;
   for i=1:Q_num
       S_xyz_temp=Q_array(i+1,1:3)';
       [S_Aoa,index_Rd]=Select_Rd(S_xyz_temp);                       % 获取方位角/俯仰角 
       StructDate{1,7}(i+1,7:8)=S_Aoa';                              % 为结构体更新对应散射体的方位角/俯仰角
       TDOA_FDOA_Az_DOA{1,8}(i+1,7:8)=S_Aoa';  
       N=20;                                                         % CFAR中每次提取当前数据点及左右长度10
       k_cfar=floor(N*3/4);
       Protect_unit = 2;                                             % Rd_数据必须是行序列
       scampling =Caiyang*1e6;
       if scampling ==Scample_rate_1
           N_pf=200;
       elseif scampling ==Scample_rate_2
           N_pf=150;
       elseif scampling ==Scample_rate_3
           N_pf=70;
       elseif scampling ==Scample_rate_4
           N_pf=30;
       elseif scampling ==Scample_rate_5
       end
       % 取出散射体对应扇区数据
       [Rd_cell{1,i},L_k]=AcquireRdDate2(RD_data,index_Rd,flag_caiyangl);  % 开始提取RD数据容器
            
       if flag_caiyang==1   % 高采样率下的分支
         Rd_array=Rd_cell{1,i};
         feng = zeros(1,L_k(2));
           for j=1:2*L_k(2)+1
               feng(j)=sum(abs(Rd_array(:,j)));   % 选取峰值所在列
           end
           [~,cfar_index]=max(feng);
           flag_noise=NoisType(abs(Rd_array(:,cfar_index)'),0.0001);      % 获取噪声类型  水平下降因子内置为0.1
       
           [Order_l,T_Sn,Xk]=Qin_CFAR(Pf,abs(Rd_array(:,cfar_index)'),2*L_k(1)+1,flag_noise,N_pf,N,k_cfar,Protect_unit);   
                                                                          % 参数类型 Pf，Rd_单列/单行，RD数据长度 噪声类型 窗口长度
                                                                          % k=（3/4-1/2）N 保护单元长度
           record_cfarout{i} = Order_l;
           record_noisetype{i}=flag_noise;
           
           NoisePower=[Xk,T_Sn(2)];
           [pufeng_array,Index_array]=RD_pufengsearch(abs(Rd_array),Order_l, NoisePower, doa_ref,i,flag_caiyangl, flag_noise); 
                                                                          %  RD矩阵  挑选出的序号
                                                                          %  噪声功率(Sn可用) 方向角 散射体序号
                                                                          %  采样率具体  噪声类型
           if pufeng_array~=-1                                            %  检查pufeng_array的有效性
               Pefeng_cell{ count_avail+1}=pufeng_array;                          
               Index_array(7,:)=index_Rd;                                 % Index_array 是列数据用于保存通过筛选的点迹信息
               TDOA_FDOA_Az_DOA{1,11+i}=Index_array';
               count_avail=count_avail+1;
           else
               Q=Q-1;
               TDOA_FDOA_Az_DOA{1,11+i}=0;                                % 填装当前时刻 RD谱检测结果
           end        
       else                                                               % 低采样率下的分支
           Rd_array=Rd_cell{1,i};
           feng = zeros(1,L_k(1));
           for j=1:2*L_k(1)+1
               feng(j)=sum(abs(Rd_array(j,:)));                           % 选取峰值所在行
           end
           [~,cfar_index]=max(feng);
           flag_noise=NoisType(abs(Rd_array(cfar_index,:)),0.0001);       % 获取噪声类型  水平下降因子内置为0.1
           [Order_l,T_Sn,Xk]=Qin_CFAR(Pf,abs(Rd_array(cfar_index,:)),2*L_k(2)+1,flag_noise,N_pf,N,k_cfar,Protect_unit);  %%参数类型 Pf，Rd_单列/单行，RD数据长度 噪声类型 窗口长度 k=（3/4-1/2）N 保护单元长度
           
           record_cfarout{i} = Order_l;
           record_noisetype{i}=flag_noise;
           
           NoisePower=[Xk,T_Sn(2)];
           [pufeng_array,Index_array]=RD_pufengsearch(abs(Rd_array),Order_l,NoisePower,doa_ref,i,flag_caiyangl,flag_noise); %%RD矩阵 挑选出的序号 噪声功率(Sn可用) 方向角 散射体序号 采样率具体
           if pufeng_array~=-1
               Pefeng_cell{count_avail+1}=pufeng_array;                          
               Index_array(7,:)=index_Rd;                                % Index_array 是列数据用于保存通过筛选的点迹信息
               TDOA_FDOA_Az_DOA{1,11+i}=Index_array';
               count_avail=count_avail+1;
           else
               Q=Q-1;
               TDOA_FDOA_Az_DOA{1,11+i}=0;
           end
       end

   end 
      format1=('TDOA-FDOA-Azimuth-DOA-%.3fMHZ-%.3f采样率-%d年%d月%d日-%d时%d分%d秒.mat');
      str1 = ['E:\混合数据帧2解析\TDOA_FDOA_Azimuth_DOA\',sprintf(format1,StructDate{2}(1),StructDate{2}(3),StructDate{1}(1),StructDate{1}(2),StructDate{1}(3),StructDate{1}(4),StructDate{1}(5),StructDate{1}(6))];
      save(str1,'TDOA_FDOA_Az_DOA')  % 保存该结构体
      str2 = sprintf('  情景2: 搜索%d个散射体，%d散射体扇区存在谱峰搜索结果！',Q_num,Q);
%     disp(['情景2状态Mention: ',sprintf('%d个升空散射体',Q) ,'谱峰搜索结果具有效数据']);
      UI_printInf{1} = record_cfarout;
      UI_printInf{2} = str2;   
      UI_printInf{3} = str1;                                              % 记录TDOA_FDOA_Az_D0A保存位置
       
    % 第三部分 参数开始遍历配对 TDOA/DOA FDOA/DOA混合定位阶段 %
    f_array=ones(1,Q_num)*StructDate{1,2}(1,1);  
    % 定位结果有两种情况： -1：无有效扇区数据   0：无适合定位点 (此时有两种可能 谱峰结果无过门限组合 或者 LMS无定位结果)
    if Q>0
       point_cluster_cell=trace_couple_index(Pefeng_cell,Q,doa_ref);    % 返回参数配对结果 元胞矩阵 矩阵行数据:
                                                                        % TDOA时差、等效信噪比_tdoa、FDOA时差、等效信噪比_fdoa、方向角_ref，俯仰角_ref 散射体索引
       if flag_caiyang==1                                               % 中高采样率定位
         if flag_caiyangl~=3
          point_cluster_cell2 = pro_trace_couple(point_cluster_cell,T_Sn(1),flag_caiyang);
         else
          point_cluster_cell2 = pro_trace_couple(point_cluster_cell,T_Sn(2),flag_caiyang);   
         end
         [loc_shunshi,couple_array] = LMS_traceloc(point_cluster_cell2,S_xyz,S_vxyz,f_array,doa_ref,flag_caiyang,flag_yundong); 
                                                                       % point_cluster,S_xyz,S_varray,f_array,doa_ref,flag_caiyang,flag_yundong
       else    %% 低采样率定位
         point_cluster_cell2 = pro_trace_couple(point_cluster_cell,T_Sn(2),flag_caiyang); 
         [loc_shunshi,couple_array] = LMS_traceloc(point_cluster_cell2,S_xyz,S_vxyz,f_array,doa_ref,flag_caiyang,flag_yundong);
       end
       
       if loc_shunshi(1,1)~=-1  %定位有结果： 有两种情况： -1是无数据 / 有效数据   couple_array 有效的配对定位组合 
           StructDate{1,11}=loc_shunshi;  
%          disp(['情景2状态Mention: ',sprintf('%d个升空散射体',Q) ,'获得有效的定位结果']);
           flag_isUsed=1;        % 数据结果可用
       else
           StructDate{1,11}=0;  
%            disp(['情景2状态Mention: ',sprintf('%d个升空散射体',Q) ,'获得无效的定位结果']);
%            disp(['异常ention: 检查point_cluster_cell2是否全0 否则 定位过程']);
           flag_isUsed=0;        % 数据结果不可用
       end
     UI_printInf{4} = couple_array;  
    else
        Loc_dingwei=-1;
        StructDate{1,11}=Loc_dingwei;
        flag_isUsed=0;           % 数据结果不可用
    end
     format1=('InstantPosition-%.3fMHZ-%.3f采样率-%d年%d月%d日-%d时%d分%d秒.mat');
     str2 = ['E:\混合数据帧2解析\CPI时刻的瞬时定位结果\定位情景2\',sprintf(format1,StructDate{2}(1),StructDate{2}(3),StructDate{1}(1),StructDate{1}(2),StructDate{1}(3),StructDate{1}(4),StructDate{1}(5),StructDate{1}(6))];
     save(str2,'StructDate')  % 保存该结构体
     UI_printInf{5} = str2;  
end  
function [Rd,L_k]=AcquireRdDate2(RD_data,index_Rd,flag_caiyangl)
    if flag_caiyangl==1
        l=163;k=1;
    elseif flag_caiyangl==2
        l=82; k=2;
    elseif flag_caiyangl==3
        l=41; k=3;
    elseif flag_caiyangl==4
        l=21; k=6;
    elseif flag_caiyangl==0
        l=1;  k=113;
    end
    L_k=[l,k];
    Rd=RD_data{index_Rd};
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 仿真必要函数 %%%%%%%%%%%%%%%%%%%%%%%%
% function StructDate = Simu_structdata()
% StructDate=cell(1,11);
% StructDate{1}=[2020,10,2,11,32,40];
% StructDate{2}=[15,30,32.5];
% StructDate{3}=[1,163,30.77];
% StructDate{4}=[0,0,0];
% StructDate{5}=[0,1];
% StructDate{6}=[0 0 0 0 0 0       %UCA坐标模块 经度 纬度 高度 速度x y z
%                0 0 0 0 0 0];     %磁x 磁y 磁z 俯仰 横向 0
% StructDate{7}=[ 3  0  0 0 0 0  0  0
%                 300,70,50 0 0 0 13.1 9.21     %% 未测向散射体Q (NumQ+1)*8  第一行：数量
%                -400 70 60 0 0 0 170.1 8.40    %% S经度 S纬度 S高度  S_速度x S_速度y S_速度z 方向 俯仰
%                 5 -100 50 0 0 0 272.9 26.53];
% StructDate{8}= [4 0 0                   %% DOA信息 (Num_DOA+1)*3 第一行 DOA的数量
%                 45 4 0 
%                 13.1 9.21 0
%                 170.1 8.40 0
%                 272.9 26.53 0           %% 方向角 俯仰角 来波方向强度                         
%                 ];
% StructDate{9}= [3 0 0 0 0
%                 5 -100 50 3 2
%                 300 70 50 1 3
%                -400 70 60 2 4];         %% zeros 1*5  T散射体信息 (Num_T+1)*5
% StructDate{10}=[45 4 1];
% end


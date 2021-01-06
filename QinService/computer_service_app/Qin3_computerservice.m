function [StructDate,UI_printInf,flag_isUsed]=Qin3_computerservice(Sigframe_decoded,Pf,flag_yundong)
%%% 情景3 Q>0 T>0
%%% Sigframe_decoder :计时器收集到的有效帧数据解析后的结果 默认是一帧元胞数据
%%% flag_yundong: 低采样率针对的是运动目标 0静止 1：开启
%%% Q、T数量 用于选择分支
%%% 输出参数： 
%%%  Loc_dingwei:    输出定位结果
%%%  StructDate:     定位结构体 某分支下
%%%  flag_isUsed=0   定位结果是否有效

%%% 特征滤波 期望信号再确认
    %%% there? or other?
addpath 'E:\A_Matlab2020a\Matlab2020a\bin\computer_service_app\ComSer2_situtation3'

% 截胡准备区
% simulation_sigframe = load('E:\混合数据帧2解析\Mixframe2_2020_11_12_10_53_58_650_decode.mat');
% simulation_sigframe =simulation_sigframe.Struct_data;
% ServiceFram =load('E:\混合数据帧2解析\Service2_Qin3Frame');
% simulation_RD=ServiceFram.RD_cell;
%

% 第一部分 数据解帧 %
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
index =1;
for i=1:P_num
    StructDate{8}(i+1,1:3) = [Sigframe_decoded{41+index},Sigframe_decoded{41+index+1},0];
    index = index+2;
end
StructDate{9}=zeros(5,5);                   % zeros 5*5  T散射体信息
T_num = Sigframe_decoded{56};
StructDate{9}(1,1) = T_num;

refdoa_order = Sigframe_decoded{65};
refS_order =1;
index=56;
for i=1:T_num
    S_order =Sigframe_decoded{index+1};
    doa_order =Sigframe_decoded{index+2};
    if doa_order == refdoa_order
        refS_order = S_order;              % 获得参考散射体的序号
    end
    StructDate{9}(i+1,1:3)=StructDate{7}(S_order+1,1:3);
    StructDate{9}(i+1,4)=S_order;
    StructDate{9}(i+1,5)=doa_order;
    index =index+2;
end

StructDate{10}=zeros(1,3);                  % 参考通道信号来波方向 DOA信息中的字段
StructDate{10}(1:2) = StructDate{8}(refdoa_order+1,1:2);
StructDate{10}(3) =  refdoa_order;
doa_ref=StructDate{1,10}(1:2);               % 参考通道DOA
%%%

if(T_num>1)
   flag_branch3=1;                                     %% 启动分支2的另一个标志
else
   flag_branch3=2;
end
TDOA_FDOA_Az_DOA_1=StructDate;                         %%TDOA-DOA-Az-FDOA数据结构体 AZ参考通道的序号
TDOA_FDOA_Az_DOA_2=StructDate;
TDOA_FDOA_Az_DOA_3=StructDate;
%%%
Q=Q_num-1;                                             % 记录搜索的扇区
T=T_num-1;                                             % 记录搜索的扇区
%%% 获取直达波方向 读取升空散射体方位扇区RD数据块
    flag_yundong=0;
    Scample_rate_1 = 32500000;        % scamping rate_1: 32.5MHZ
    Scample_rate_2 = 16250000;        % scamping rate_2: 16.25MHZ
    Scample_rate_3 = 8125000;         % scamping rate_3: 8.125MHZ
    Scample_rate_4 = 4062500;         % scamping rate_4: 4.0625MHZ
    Scample_rate_5 = 203125;          % scamping rate_5: 203125KHZ       
%%%

%%%% 存储TDOA-DOA-Azimuth-DOA结构体
    for i=10:-1:5
        TDOA_FDOA_Az_DOA_1{i+1}=TDOA_FDOA_Az_DOA_1{i};             %% 空出胞单元5 填充方位扇区分布信息
        TDOA_FDOA_Az_DOA_2{i+1}=TDOA_FDOA_Az_DOA_2{i};             %% 空出胞单元5 填充方位扇区分布信息
        TDOA_FDOA_Az_DOA_3{i+1}=TDOA_FDOA_Az_DOA_3{i};             %% 空出胞单元5 填充方位扇区分布信息
    end
        TDOA_FDOA_Az_DOA_1{5}=[1,12,30];                           %% 方位扇区的序号范围、方位分辨率
        TDOA_FDOA_Az_DOA_2{5}=[1,12,30];                           %% 方位扇区的序号范围、方位分辨率
        TDOA_FDOA_Az_DOA_3{5}=[1,12,30];                           %% 方位扇区的序号范围、方位分辨率
%%%

%%% 获取IQ信号采样率和带宽
    Caiyang=StructDate{1,2}(1,3);
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
   Q_array=StructDate{1,7};  %%% 提取出散射体的数量 经/纬/高 瞬时速度
   T_array=StructDate{1,9};  %%% 5*5 散射体位置信息 在升空散射体字段中的序号 在DOA估计信息字段中的序号
   
   %%% 获取散射体坐标、获取散射体速度
    cen_state=StructDate{6}(1,1:3); %% 获取UCA经纬高坐标
  % true area
    S_xyz = LBM_sXYZ((Q_array(2:Q_num+1,1:3))',cen_state);
    S_xyz = Q_array(2:Q_num+1,1:3)';               %  自制的解析帧散射体数据是站心坐标系 不需要坐标转换         
  % temple area
   % S_xyz = StructDate{1,7}(2:end,1:3)';
    

%%%     true area                 % 截胡区域关闭则将此区域打开
    S_ref =S_xyz(:,refS_order);   % 默认将第一个可观测的散射体作为参考散射体  修改:默认将最强散射体的来波方向作为参考 
    T_sxyz=zeros(3,5);            
    j=1;
    for i=1:T_num                 % 获取不包含参考散射体的所有散射体坐标
        if T_array(i+1,4) ~= refS_order
            T_sxyz(:,j)=S_xyz(:,T_array(i+1,4));
            j=j+1;
        end
    end

    vsi_array=Q_array(2:Q_num+1,4);
    fai_array=Q_array(2:Q_num+1,7);
    UCA_sudu =[0 0];                                                 %% UCA速度未知     
    S_vxyz=Count_Sv(S_xyz,vsi_array,fai_array,UCA_sudu);             %% 计算散射体数据 散射体坐标、散射体瞬时速度、航向角、UCA速度1*2(暂无)
    S_vref=S_vxyz(:,T_array(2,4));                                   %% 记录下参考散射体的速度
    
    StructDate{1,7}(2:Q_num+1,4:6)=S_vxyz';                          %% 散射体速度赋给结构体
    TDOA_FDOA_Az_DOA_1{1,8}(2:Q_num+1,4:6)=S_vxyz';
    TDOA_FDOA_Az_DOA_2{1,8}(2:Q_num+1,4:6)=S_vxyz';
    TDOA_FDOA_Az_DOA_3{1,8}(2:Q_num+1,4:6)=S_vxyz';    
    RD_data =cell(1,12);                                            % Sigframe_decoded 提取RD谱数据
    for i=68:79
        RD_data{i-67} = Sigframe_decoded{i};
    end
    
    % 截胡 %
%     StructDate=Simu_structdata();
%      RD_data = simulation_RD;
%     S_xyz=StructDate{1,7}(2:Q_num+1,1:3)';
%     S_vxyz=Count_Sv(S_xyz,vsi_array,fai_array,UCA_sudu);             %% 计算散射体数据 散射体坐标、散射体瞬时速度、航向角、UCA速度1*2(暂无)
%     StructDate{1,7}(2:Q_num+1,4:6)=S_vxyz';
%     doa_ref=StructDate{1,10}(1:2);                                   % 参考通道DOA
%     refS_order = 3;
%     S_ref =S_xyz(:,refS_order);   %% 默认将第一个可观测的散射体作为参考散射体  修改:默认将最强散射体的来波方向作为参考 
%     T_array=StructDate{1,9};  %%% 5*5 散射体位置信息 在升空散射体字段中的序号 在DOA估计信息字段中的序号
%     T_sxyz=zeros(3,5); 
%     j=1;
%     for i=1:T_num                 % 获取不包含参考散射体的所有散射体坐标
%         if T_array(i+1,4) ~= refS_order
%             T_sxyz(:,j)=S_xyz(:,T_array(i+1,4));
%             j=j+1;
%         end
%     end
    %
    UI_printInf =cell(1,11);
%%%

%%% 第二部分 获取散射体数量 散射体对应方位角/俯仰角 填装Rd_cell CFAR检测 峰值检测   
   Rd1_cell=cell(1,6);       %%% 声明Rd1数据容器  用于装载分支1 p-1个DOA对应的方位扇区的RD数据
   Rd2_cell=cell(1,3);       %%% 声明Rd2数据容器  用于装载分支2 Q-1个非参考散射体对应方位扇区RD数据
   Rd3_cell=cell(1,3);       %%% 声明Rd3数据容器  用于装载分支3 T-1个非参考散射体对应方位扇区RD数据
   Pefeng1_cell=cell(1,6);   %%% Pefeng1_cell容器 用于装载分支1 p-1个DOA有效的谱峰搜索结果
   Pefeng2_cell=cell(1,3);   %%% Pefeng2_cell容器 用于装载分支2 Q-1个散射体有效的谱峰搜索结果
   Pefeng3_cell=cell(1,3);   %%% Pefeng3_cell容器 用于装载分支2 T-1个可观测有效的谱峰搜索结果

    
   %%% CFAR 检测器需要的参数
   N=20;
   k=floor(N*3/4);
   Protect_unit = 2;  %% Rd_数据必须是行序列
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
  
   % 根据分支生成峰值搜索结果   
      Order_Sdoaref=refdoa_order;          %% 参考T的DOA序号    必须保证T_array 第一行是参考散射体
      DOA_array=StructDate{1,8};
      DOA_num = DOA_array(1,1)-1;          %% 记录除去参考散射体方向以外的DOA角度
      
       T_doa =zeros(2,5);  T_v=zeros(3,5);
       j=1;
       for i=1:T_num                       % 获取不包含参考散射体的方向角和速度
        if T_array(i+1,4)~=refS_order
        T_doa(:,j)=(DOA_array(T_array(i+1,5)+1,1:2))';                  %% DOA
        T_v(:,j)=S_vxyz(:,T_array(i+1,4));
        j=j+1;
        end
       end
        S_vxyz(:,refS_order)=[];                                      % 从中删除参考散射体的速度
        S_xyz(:,refS_order) = [];                                     % 从中删除参考散射体
        DOA_array(Order_Sdoaref+1,:)=[];                              % 删除参考散射体的DOA行 第一行是DOA的数量
        DOA_num2 = DOA_num-1;                                           % 用于记录有效的DOA扇区
        doa_array =DOA_array(2:end,1:2)';

      % 开始提取p-1个RD数据谱 
      rd_index1=zeros(1,DOA_num);  
      count_branch1=0;
      record_cfarout_branch1 = cell(1,DOA_num);
      for i=1:DOA_num                              % 分支1
         rd_index1(i) =Select_RdBaseDoa(DOA_array(i+1,1));
         [Rd1_cell{1,i},L_k]=AcquireRdDate2(RD_data,rd_index1(i),flag_caiyangl); % 开始提取RD数据容器
          
         if flag_caiyang==1                        % 高采样率下的分支
           Rd_array=Rd1_cell{1,i};
           feng = zeros(1,L_k(2));
           for j=1:2*L_k(2)+1
               feng(j)=sum(abs(Rd_array(:,j)));    % 选取峰值所在列
           end
           [~,cfar_index]=max(feng);
           flag_noise=NoisType(abs(Rd_array(:,cfar_index)'),0.0001); % 获取噪声类型 0.0001是置信水平 水平下降因子内置为0.1
           [Order_l,T_Sn,Xk]=Qin_CFAR(Pf,abs(Rd_array(:,cfar_index)'),2*L_k(1)+1,flag_noise,N_pf,N,k,Protect_unit);   
                                                                     % 参数类型 Pf，Rd_单列/单行，RD数据长度 噪声类型 窗口长度
                                                                     % k=（3/4-1/2）N 保护单元长度   
           record_cfarout_branch1{i} = Order_l;
           NoisePower=[Xk,T_Sn(2)];
           [pufeng_array,Index_array]=RD_pufengsearch(abs(Rd_array), Order_l,   NoisePower,   doa_ref,  i,   flag_caiyangl, flag_noise); 
                                                                     %  RD矩阵  挑选出的序号
                                                                     %  噪声功率(Sn可用) 方向角 散射体序号
                                                                     %  采样率具体  噪声类型
           if pufeng_array~=-1                                       % 检查pufeng_array的有效性
               Pefeng1_cell{count_branch1+1}=pufeng_array;                          
               Index_array(7,:)=rd_index1(i);                        % Index_array 是列数据用于保存通过筛选的扇区
            %  doa_array(:,count_branch1+1)=(DOA_array(i+1,1:2))';
               TDOA_FDOA_Az_DOA_1{1,11+i}=Index_array';
               count_branch1=count_branch1+1;
           else
               DOA_num2=DOA_num2-1;
               TDOA_FDOA_Az_DOA_1{1,11+i}=0;                           %% 填装当前时刻 RD谱检测结果
           end
           
         else                                   % 低采样率下的分支
           Rd_array=Rd1_cell{1,i};
           feng = zeros(1,L_k(1));
           for j=1:2*L_k(1)+1
               feng(j)=sum(abs(Rd_array(j,:)));  % 选取峰值所在行
           end
           [~,cfar_index]=max(feng);
           flag_noise=NoisType(abs(Rd_array(cfar_index,:)),0.0001);    %% 获取噪声类型  水平下降因子内置为0.1
           [Order_l,T_Sn,Xk]=Qin_CFAR(Pf,abs(Rd_array(cfar_index,:)),2*L_k(2)+1,flag_noise,N_pf,N,k,Protect_unit);  %%参数类型 Pf，Rd_单列/单行，RD数据长度 噪声类型 窗口长度 k=（3/4-1/2）N 保护单元长度
           record_cfarout_branch1{i} = Order_l;
           
           NoisePower=[Xk,T_Sn(2)];
           [pufeng_array,Index_array]=RD_pufengsearch(abs(Rd_array),Order_l,NoisePower,doa_ref,i,flag_caiyangl,flag_noise); %%RD矩阵 挑选出的序号 噪声功率(Sn可用) 方向角 散射体序号 采样率具体
           if pufeng_array~=-1
               Pefeng1_cell{count_branch1+1}=pufeng_array;                          
               Index_array(7,:)=rd_index1(i);                            %% Index_array 是列数据用于保存通过筛选的点迹信息
             % doa_array(:,count_branch1+1)=(DOA_array(i+1,1:2))';
               TDOA_FDOA_Az_DOA_1{1,11+i}=Index_array';
               count_branch1=count_branch1+1;
           else
               DOA_num2=DOA_num2-1;
               TDOA_FDOA_Az_DOA_1{1,11+i}=0;
           end
         end   
      end
      format1=('TDOA-FDOA-Azimuth-DOA-%.3fMHZ-%.3f采样率-%d年%d月%d日-%d时%d分%d秒.mat');
      str1 = ['E:\混合数据帧2解析\TDOA_FDOA_Azimuth_DOA_sit3_1\',sprintf(format1,StructDate{2}(1),StructDate{2}(3),StructDate{1}(1),StructDate{1}(2),StructDate{1}(3),StructDate{1}(4),StructDate{1}(5),StructDate{1}(6))];
      save(str1,'TDOA_FDOA_Az_DOA_1')  %%保存该结构体
%       disp(['情景3分支1: ',sprintf('共搜索%d个扇区，%d组DOA值对应扇区',DOA_num-1,DOA_num2) ,'存在谱峰搜索结果']);
      str2 = sprintf('  情景3分支1: 搜索%d个扇区，%d扇区存在谱峰结果！',DOA_num-1,DOA_num2);
      UI_printInf{1} =record_cfarout_branch1;
      UI_printInf{2} = str2;
      UI_printInf{3} = str1;
            
   if ((Q_num>1 && T_num==1) || flag_branch3==1) %% 分支2
     count_branch2=0;
     record_cfarout_branch2 = cell(1,Q_num-1);                                  % 记录CFAR输出结果送入到UI中
     
     for i=1:(Q_num-1)
       [S_Aoa,index_Rd]=Select_Rd(S_xyz(:,i));                                % 获取方位角/俯仰角    
       StructDate{1,7}(i+1,7:8)=S_Aoa';                                       % 为结构体更新对应散射体的方位角/俯仰角
       TDOA_FDOA_Az_DOA_2{1,8}(i+1,7:8)=S_Aoa';
       [Rd2_cell{1,i},L_k]=AcquireRdDate2(RD_data,index_Rd,flag_caiyangl);    % 开始提取RD数据容器
       
       if flag_caiyang==1   %%% 高采样率下的分支
           Rd_array=Rd2_cell{1,i};
           feng = zeros(1,L_k(2));
           for j=1:2*L_k(2)+1
               feng(j)=sum(abs(Rd_array(:,j)));   %% 选取峰值所在列
           end
           [~,cfar_index2]=max(feng);
           flag_noise=NoisType(abs(Rd_array(:,cfar_index2)'),0.0001);                   %% 获取噪声类型  0.0001是置信水平 水平下降因子内置为0.1
           [Order_l,T_Sn,Xk]=Qin_CFAR(Pf,abs(Rd_array(:,cfar_index2)'),2*L_k(1)+1,flag_noise,N_pf,N,k,Protect_unit);   
                                           %%% 参数类型 Pf，Rd_单列/单行，RD数据长度 噪声类型 窗口长度
                                           %%% k=（3/4-1/2）N 保护单元长度  
           record_cfarout_branch2{i} = Order_l;
           
           NoisePower=[Xk,T_Sn(2)];
           [pufeng_array,Index_array]=RD_pufengsearch(abs(Rd_array), Order_l,   NoisePower,   doa_ref,  i,   flag_caiyangl, flag_noise); 
                                                                     %%%  RD矩阵  挑选出的序号
                                                                     %%%  噪声功率(Sn可用) 方向角 散射体序号
                                                                     %%%  采样率具体  噪声类型
           if pufeng_array~=-1             %%  检查pufeng_array的有效性
               Pefeng2_cell{count_branch2+1}=pufeng_array;                          
               Index_array(7,:)=index_Rd;                            %% Index_array 是列数据用于保存通过筛选的点迹信息
               TDOA_FDOA_Az_DOA_2{1,11+i}=Index_array';
               count_branch2=count_branch2+1;
           else
               Q=Q-1;
               TDOA_FDOA_Az_DOA_2{1,11+i}=0;                           %% 填装当前时刻 RD谱检测结果
           end
           
       else   % 低采样率下的分支
           Rd_array=Rd2_cell{1,i};
           feng = zeros(1,L_k(1));
           for j=1:2*L_k(1)+1
               feng(j)=sum(abs(Rd_array(j,:)));   %% 选取峰值所在行
           end
           [~,cfar_index2]=max(feng);
           flag_noise=NoisType(abs(Rd_array(cfar_index2,:)),0.0001);  % 获取噪声类型 水平下降因子内置为0.1
           [Order_l,T_Sn,Xk]=Qin_CFAR(Pf,abs(Rd_array(cfar_index2,:)),2*L_k(2)+1,flag_noise,N_pf,N,k,Protect_unit);  %%参数类型 Pf，Rd_单列/单行，RD数据长度 噪声类型 窗口长度 k=（3/4-1/2）N 保护单元长度
           record_cfarout_branch2{i} = Order_l;                        % 记录CFAR输出结果送入到UI中
           
           NoisePower=[Xk,T_Sn(2)];
           [pufeng_array,Index_array]=RD_pufengsearch(abs(Rd_array),Order_l,NoisePower,doa_ref,i,flag_caiyangl,flag_noise); %%RD矩阵 挑选出的序号 噪声功率(Sn可用) 方向角 散射体序号 采样率具体
           if pufeng_array~=-1
               Pefeng2_cell{count_branch2+1}=pufeng_array;
               Index_array(7,:)=index_Rd;                             % Index_array 是列数据
               TDOA_FDOA_Az_DOA_2{1,11+i}=Index_array';
               count_branch2=count_branch2+1;
           else
               Q=Q-1;
               TDOA_FDOA_Az_DOA_2{1,11+i}=0;
           end
       end
     end
      format1=('TDOA-FDOA-Azimuth-DOA-%.3fMHZ-%.3f采样率-%d年%d月%d日-%d时%d分%d秒.mat');
      str1 = ['E:\混合数据帧2解析\TDOA_FDOA_Azimuth_DOA_sit3_2\',sprintf(format1,StructDate{2}(1),StructDate{2}(3),StructDate{1}(1),StructDate{1}(2),StructDate{1}(3),StructDate{1}(4),StructDate{1}(5),StructDate{1}(6))];
      save(str1,'TDOA_FDOA_Az_DOA_2')  %%保存该结构体
      str2 = sprintf('  情景3分支2: 搜索%d个散射体，%d散射体扇区存在谱峰结果！',Q_num-1,Q);
%       disp(['情景3状态Mention分支2: ',sprintf('%d个散射体',Q) ,'谱峰搜索结果有效']);
      UI_printInf{4} =record_cfarout_branch2;
      UI_printInf{5} = str2;
      UI_printInf{6} = str1;

   end  
   
   if (Q_num>1 && T_num>1)  % 分支3
      count_branch3=0;
      record_cfarout_branch3 = cell(1,T_num-1);
      for i=1:(T_num-1)
       index_Rd=Select_RdBaseDoa(T_doa(1,i));                                 % 获取方位角/俯仰角 
       
       StructDate{1,7}(i+1,7:8)=S_Aoa';                                       % 为结构体更新对应散射体的方位角/俯仰角
       TDOA_FDOA_Az_DOA_3{1,8}(i+1,7:8)=T_doa(:,i)';
       [Rd3_cell{1,i},L_k]=AcquireRdDate2(RD_data,index_Rd,flag_caiyangl);    % 开始提取RD数据容器
       if flag_caiyang==1   % 高采样率下的分支
           Rd_array=Rd3_cell{1,i};
           feng = zeros(1,L_k(2));
           for j=1:2*L_k(2)+1
               feng(j)=sum(abs(Rd_array(:,j)));   % 选取峰值所在列
           end
           [~,cfar_index3]=max(feng);
           flag_noise=NoisType(abs(Rd_array(:,cfar_index3)'),0.0001);          % 获取噪声类型 水平下降因子内置为0.1
           [Order_l,T_Sn,Xk]=Qin_CFAR(Pf,abs(Rd_array(:,cfar_index3)'),2*L_k(1)+1,flag_noise,N_pf,N,k,Protect_unit);   
                                           %%% 参数类型 Pf，Rd_单列/单行，RD数据长度 噪声类型 窗口长度
                                           %%% k=（3/4-1/2）N 保护单元长度
           record_cfarout_branch3{i} = Order_l;
           
           NoisePower=[Xk,T_Sn(2)];
           [pufeng_array,Index_array]=RD_pufengsearch(abs(Rd_array), Order_l,   NoisePower,   doa_ref,  i,   flag_caiyangl, flag_noise); 
                                                                     %%%  RD矩阵  挑选出的序号
                                                                     %%%  噪声功率(Sn可用) 方向角 散射体序号
                                                                     %%%  采样率具体  噪声类型
           if pufeng_array~=-1              %%  检查pufeng_array的有效性
               Pefeng3_cell{count_branch3+1}=pufeng_array;                          
               Index_array(7,:)=index_Rd;                            %% Index_array 是列数据用于保存通过筛选的点迹信息
               TDOA_FDOA_Az_DOA_3{1,11+i}=Index_array';
               count_branch3=count_branch3+1;
           else
               T=T-1;
               TDOA_FDOA_Az_DOA_3{1,11+i}=0;                        % 填装当前时刻 RD谱检测结果
           end
       else                 %%% 低采样率下的分支
           Rd_array=Rd3_cell{1,i};
           feng = zeros(1,L_k(1));
           for j=1:2*L_k(1)+1
               feng(j)=sum(abs(Rd_array(j,:)));   % 选取峰值所在行
           end
           [~,cfar_index3]=max(feng);
           flag_noise=NoisType(abs(Rd_array(cfar_index3,:)),0.0001); % 获取噪声类型  水平下降因子内置为0.1
           [Order_l,T_Sn,Xk]=Qin_CFAR(Pf,abs(Rd_array(cfar_index3,:)),2*L_k(2)+1,flag_noise,N_pf,N,k,Protect_unit);  %%参数类型 Pf，Rd_单列/单行，RD数据长度 噪声类型 窗口长度 k=（3/4-1/2）N 保护单元长度
           record_cfarout_branch3{i} = Order_l;
           
           NoisePower=[Xk,T_Sn(2)];
           [pufeng_array,Index_array]=RD_pufengsearch(abs(Rd_array),Order_l,NoisePower,doa_ref,i,flag_caiyangl,flag_noise); %%RD矩阵 挑选出的序号 噪声功率(Sn可用) 方向角 散射体序号 采样率具体
           if pufeng_array~=-1
               Pefeng3_cell{count_branch3+1}=pufeng_array;
               Index_array(7,:)=index_Rd;                            % Index_array 是列数据
               TDOA_FDOA_Az_DOA_3{1,11+i}=Index_array';              % CFAR检测结果序号包含在其中
               count_branch3=count_branch3+1;
           else
               T=T-1;
               TDOA_FDOA_Az_DOA_3{1,11+i}=0;
           end
       end
     end
%       disp(['情景3状态Mention分支3: ',sprintf('%d个可观测散射体',T) ,'谱峰搜索结果有效']);
      format1=('TDOA-FDOA-Azimuth-DOA-%.3fMHZ-%.3f采样率-%d年%d月%d日-%d时%d分%d秒.mat');
      str1 = ['E:\混合数据帧2解析\TDOA_FDOA_Azimuth_DOA_sit3_3\',sprintf(format1,StructDate{2}(1),StructDate{2}(3),StructDate{1}(1),StructDate{1}(2),StructDate{1}(3),StructDate{1}(4),StructDate{1}(5),StructDate{1}(6))];
      save(str1,'TDOA_FDOA_Az_DOA_3')  %%保存该结构体
      str2 = sprintf('  情景3分支3: 搜索%d个散射体，%d散射体所在扇区存在谱峰结果！',T_num-1,T);
      UI_printInf{7} =record_cfarout_branch3;
      UI_printInf{8} = str2;
      UI_printInf{9} = str1;
   end
   

   
%%%%% 第三部分 参数开始遍历配对 TDOA/DOA FDOA/DOA混合定位阶段%%%%%
    if (Q_num==1 && T_num==1)
       [trace_cell_1,is_trace1]=Q3_tracecouple_branch1(Pefeng1_cell,DOA_num2,doa_array);   %% 有效的DOA_num2不足 point_cell_1=0
       flag_branch=1;
    elseif (Q_num>1 && T_num==1)
       [trace_cell_1,is_trace1]=Q3_tracecouple_branch1(Pefeng1_cell,DOA_num2,doa_array);   %% 有效的DOA_num2不足 point_cell_1=0
       [trace_cell_2,is_trace2]=Q3_tracecouple_branch2(Pefeng2_cell,Q_num);                %% 有效的Num_Q不足    point_cell_2=0
       flag_branch=2;
    elseif (Q_num>1 && T_num>1) 
       [trace_cell_1,is_trace1]=Q3_tracecouple_branch1(Pefeng1_cell,DOA_num2,doa_array);   %% 有效的DOA_num2不足 point_cell_1=0
       [trace_cell_2,is_trace2]=Q3_tracecouple_branch2(Pefeng2_cell,Q_num);                %% 有效的Num_Q不足    point_cell_2=0
       [trace_cell_3,is_trace3]=Q3_tracecouple_branch3(Pefeng3_cell,T_num);               %% 有效的Num_T不足    point_cell_3=0
       flag_branch=3;
    end
   
    trace_array_1=0;       trace_array_2=0;               trace_array_3=0;
    if (Q_num==1 && T_num==1)
       if is_trace1~=0
           if flag_caiyangl~=3
               trace_cell2_1=Qin3_pro_trace_couple(trace_cell_1,T_Sn(1),flag_caiyang);
               trace_array_1=trace_cell2_1{1};
           else
               trace_cell2_1=Qin3_pro_trace_couple(trace_cell_1,T_Sn(2),flag_caiyang);
               trace_array_1=trace_cell2_1{1};
           end
       else
%            disp('情景3状态Mention分支1: 任务1无有效定位参数组合');
       end
       
    elseif (Q_num>1 && T_num==1)
       if is_trace1~=0
           if flag_caiyangl~=3
               trace_cell2_1=Qin3_pro_trace_couple(trace_cell_1,T_Sn(1),flag_caiyang);
               trace_array_1=trace_cell2_1{1};
           else
               trace_cell2_1=Qin3_pro_trace_couple(trace_cell_1,T_Sn(2),flag_caiyang);
               trace_array_1=trace_cell2_1{1};
           end
       else
%             disp('情景3状态Mention分支2:  任务1无有效定位参数组合');
       end
       if is_trace2~=0
           if flag_caiyangl~=3
               trace_cell2_2=Qin3_pro_trace_couple(trace_cell_2,T_Sn(1),flag_caiyang);
               trace_array_2=trace_cell2_2{1};
           else
               trace_cell2_2=Qin3_pro_trace_couple(trace_cell_2,T_Sn(2),flag_caiyang);
               trace_array_2=trace_cell2_2{1};
           end
       else
%              disp('情景3状态Mention分支2:  任务2无有效定位参数组合');
       end
       
    elseif (Q_num>1 && T_num>1)
       if is_trace1~=0
           if flag_caiyangl~=3
               trace_cell2_1=Qin3_pro_trace_couple(trace_cell_1,T_Sn(1),flag_caiyang);
               trace_array_1=trace_cell2_1{1};
           else
               trace_cell2_1=Qin3_pro_trace_couple(trace_cell_1,T_Sn(2),flag_caiyang);
               trace_array_1=trace_cell2_1{1};
           end
       else
%             disp('情景3状态Mention分支3:  任务1无有效定位参数组合');
       end
       if is_trace2~=0
           if flag_caiyangl~=3
               trace_cell2_2=Qin3_pro_trace_couple(trace_cell_2,T_Sn(1),flag_caiyang);
               trace_array_2=trace_cell2_2{1};
           else
               trace_cell2_2=Qin3_pro_trace_couple(trace_cell_2,T_Sn(2),flag_caiyang);
               trace_array_2=trace_cell2_2{1};
           end
       else
%             disp('情景3状态Mention分支3:  任务2无有效定位参数组合');
       end
       if is_trace3~=0
           if flag_caiyangl~=3
               trace_cell2_3=Qin3_pro_trace_couple(trace_cell_3,T_Sn(1),flag_caiyang);
               trace_array_3=trace_cell2_3{1};
               
           else
               trace_cell2_3=Qin3_pro_trace_couple(trace_cell_3,T_Sn(2),flag_caiyang);
               trace_array_3=trace_cell2_3{1};
           end
       else
%             disp('情景3状态Mention分支3: 任务3无有效定位参数组合');
       end
    end
    
    trace_cluster=cell(1,3);
    trace_cluster{1}=trace_array_1; trace_cluster{2}=trace_array_2;  trace_cluster{3}=trace_array_3;
    S_txyz=[S_ref ,T_sxyz];
    S_tv  =[S_vref,T_v];
    f_array=ones(1,4)*StructDate{1,2}(1,1);
    [loc_shunshi,couple_array]=Qin3_LMStraceloc(trace_cluster,S_xyz,S_txyz,S_vxyz,S_tv,f_array,doa_ref,flag_caiyang,flag_branch,Q_num,T_num); 
    
    if loc_shunshi~=-1        % 定位有结果： 有两种情况： -1是无数据 / 有效数据   couple_array 有效的配对定位组合 
       StructDate{1,11}=loc_shunshi;  
%        disp(['情景3状态Mention: ',sprintf('分支%d',flag_branch) ,'定位结果有效']);
       flag_isUsed=1;        % 数据结果可用
    else
       StructDate{1,11}=0;  
%        disp(['情景3状态Mention: ',sprintf('分支%d',flag_branch) ,'定位结果无效']);
       disp('异常ention: 检查 trace_cluster是否全0 否则 定位过程');
       flag_isUsed=0;        % 数据结果不可用
    end
    format1=('InstantPosition-%.3fMHZ-%.3f采样率-%d年%d月%d日-%d时%d分%d秒.mat');
    str2 = ['E:\混合数据帧2解析\CPI时刻的瞬时定位结果\定位情景3\',sprintf(format1,StructDate{2}(1),StructDate{2}(3),StructDate{1}(1),StructDate{1}(2),StructDate{1}(3),StructDate{1}(4),StructDate{1}(5),StructDate{1}(6))];
    save(str2,'StructDate')  % 保存该结构体
    UI_printInf{10} = couple_array; 
    UI_printInf{11} = str2; 
end


function index_Rd=Select_RdBaseDoa(fai)
% index_Rd 是索引号
if fai<180
    if fai<90
       if fai>=0 && fai<30
           index_Rd=1;
       elseif fai>=30 && fai<60
           index_Rd=2;
       else
           index_Rd=3;
       end
    else
       if fai>=90 && fai<120
           index_Rd=4;
       elseif fai>=120 && fai<150
           index_Rd=5;
       else
           index_Rd=6;
       end
    end
else
    if fai<270
        if fai>=180 && fai<210
           index_Rd=7;
       elseif fai>=210 && fai<240
           index_Rd=8;
       else
           index_Rd=9;
       end
    else
       if fai>=270 && fai<300
           index_Rd=10;
       elseif fai>=300 && fai<330
           index_Rd=11;
       else
           index_Rd=12;
       end
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 仿真必要函数 %%%%%%%%%%%%%%%%%%%%%%%%
function [StructDate] = Simu_structdata()
StructDate=cell(1,11);
StructDate{1}=[2020,10,2,11,32,40,1];
StructDate{2}=[15,30,32.5];
StructDate{3}=[1,163,30.77];
StructDate{4}=[0,0,0];
StructDate{5}=[0,1];
StructDate{6}=[0 0 0 0 0 0       %UCA坐标模块 经度 纬度 高度 速度x y z
               0 0 0 0 0 0];     %磁x 磁y 磁z 俯仰 横向 0
StructDate{7}=[ 3  0  0 0 0 0  0  0
                300,70,50 0 0 0 13.1 9.21     %% 未测向散射体Q (NumQ+1)*8  第一行：数量
               -400 70 60 0 0 0 170.1 8.40    %% S经度 S纬度 S高度  S_速度x S_速度y S_速度z 方向 俯仰
                5 -100 50 0 0 0 272.9 26.53];
StructDate{8}= [4 0 0                   %% DOA信息 (Num_DOA+1)*3 第一行 DOA的数量
                45 4 0 
                272.9 26.53 0           %% 方向角 俯仰角 来波方向强度                         
                13.1 9.21 0
                170.1 8.40 0
                ];
StructDate{9}= [3 0 0 0 0
                5 -100 50 3 2
                300 70 50 1 3
               -400 70 60 2 4];         %% zeros 1*5  T散射体信息 (Num_T+1)*5
StructDate{10}=[272.9 26.53 2];
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
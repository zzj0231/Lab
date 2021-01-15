function [StructDate,UI_printInf,flag_isUsed]=Qin1_computerseverceFormal(SignalState_cell,Num_m,Num_zhen,ya, yb, k_up, k_down, r)
%%% 情景一计算服务器 
%  Num_m               定位站数
%  Num_zhen            站内效帧数
%  SignalState_cell    是当前站的有效帧数
%  flag_isUsed =0;     用于判断最终数据是否是有效数据
addpath(genpath('E:\A_Matlab2020a\Matlab2020a\bin\computer_service_app\ComSer2_situtation1'))

       Process_M = load('E:\混合数据帧2解析\情景1中间结果存储\Process_M.mat').Process_M;    
 ClusterDoa_cell = load('E:\混合数据帧2解析\情景1中间结果存储\ClusterDoa_cell.mat').ClusterDoa_cell;  
StateLoc_recoder = load('E:\混合数据帧2解析\情景1中间结果存储\StateLoc_recoder.mat').StateLoc_recoder;
targetLocResult  = load('E:\混合数据帧2解析\情景1中间结果存储\targetLocResult.mat').targetLocResult;        
UI_printInf = cell(1,4);
FrameInState_num = Num_zhen;
flag_isUsed = 0;
StructDate =0;

if Process_M == 0                                                                % 当是起始站 不需要利用XYZ_ref进行比较
    CPI_DOAcell = AcquireDoA(SignalState_cell,FrameInState_num);                 % CPI_DOAcell是元胞矩阵 第一帧数据的DOA已经进行过滤
    ClusterDoa_cell{Process_M+1}=DOA_cluster(CPI_DOAcell,FrameInState_num);      % 存储 功能: 收集当前站DOA聚类结果      
    StateLoc_recoder(:,Process_M+1)=[0,0,0]';                                    % 存储 功能: 记录当前站数
    LBH_ref = [SignalState_cell{1}{12},SignalState_cell{1}{13},SignalState_cell{1}{14}]';
    StateLoc_recoder(:,Process_M+1)=LBH_ref';                                    % 正式发布删除 建立以初始站位为核心的坐标系，初始站位[0,0,0]
    Process_M = Process_M+1;
    str = sprintf('  情景1：第%d站处理完毕！',Process_M);
    UI_printInf{1}=str;
    save('E:\混合数据帧2解析\情景1中间结果存储\Process_M.mat','Process_M');
    save('E:\混合数据帧2解析\情景1中间结果存储\ClusterDoa_cell.mat','ClusterDoa_cell');
    save('E:\混合数据帧2解析\情景1中间结果存储\StateLoc_recoder.mat','StateLoc_recoder');
    save('E:\混合数据帧2解析\情景1中间结果存储\LBH_ref.mat','LBH_ref');
elseif (Process_M < Num_m && Process_M>0)
    LBH_next =[SignalState_cell{1}{12},SignalState_cell{1}{13},SignalState_cell{1}{14}]';
    LBH_ref = load('E:\混合数据帧2解析\情景1中间结果存储\LBH_ref.mat').LBH_ref;
  % XYZ_next = LBH_XYZ(LBH_next,LBH_ref);   % 无真实数据时，收到的数据时站心坐标
    XYZ_next = LBH_next;
    CPI_DOAcell = AcquireDoA(SignalState_cell,FrameInState_num);        % CPI_DOAcell是元胞矩阵 第一帧数据的DOA已经进行过滤
    ClusterDoa_cell{Process_M+1} = DOA_cluster(CPI_DOAcell,FrameInState_num); % 收集当前站DOA聚类结果;
    StateLoc_recoder(:,Process_M+1) = XYZ_next;
    Process_M = Process_M+1;
    lineloc_res = line_location(StateLoc_recoder(:,Process_M-1:Process_M),ClusterDoa_cell(Process_M-1:Process_M));
    lineloc_res(:,all(lineloc_res==0,1)) = [];   
    NumResult=size(lineloc_res,2);                                              
    if NumResult>0
        targetLocResult=[targetLocResult lineloc_res];
    end
    save('E:\混合数据帧2解析\情景1中间结果存储\Process_M.mat','Process_M');
    save('E:\混合数据帧2解析\情景1中间结果存储\ClusterDoa_cell.mat','ClusterDoa_cell');
    save('E:\混合数据帧2解析\情景1中间结果存储\StateLoc_recoder.mat','StateLoc_recoder');
    save('E:\混合数据帧2解析\情景1中间结果存储\targetLocResult.mat','targetLocResult');
    disp(['情景1状态Mention：获得第',sprintf('%d',Process_M),'站有效数据']);
    str = sprintf('  情景1：第%d站处理完毕！',Process_M);
    UI_printInf{1}=str;
end
% 第二部分 减法聚类获得最终目标位置
if (Process_M >= Num_m)
    ClusterDoa_cell = cell(1,1);
    StateLoc_recoder = zeros(3,1);
    
    targetLocResult(:,all(targetLocResult==0,1))=[];
    NumResult=size(targetLocResult);
    str = sprintf("  情景1：测向交汇，获得%d组测向交汇定位结果！",NumResult(2));
    UI_printInf{2}=str;
    if NumResult(2)>0
        Loc_sub=subcluster_loc(targetLocResult,Num_m,ya, yb, k_up, k_down, r);
        dim_sub=size(Loc_sub);
        if dim_sub(2)>0
            flag_isUsed=1;
        else
            Loc_sub=0;
            flag_isUsed=0;
        end
    else
        Loc_sub=-1;
        flag_isUsed=0;
    end
    UI_printInf{4}=Loc_sub;
    loc_end=zeros(3,200);
    if  flag_isUsed==1
        for i=1:dim_sub(2)
            loc_group= Loc_sub{i};
            loc_end(1,i)=sum(loc_group(1,:))/size(loc_group,2);
            loc_end(2,i)=sum(loc_group(2,:))/size(loc_group,2);
        end
    end
    loc_end(:,all(loc_end==0,1))=[];
    str = sprintf("  情景1：获得%d组定位结果！",size(loc_end,2));
    UI_printInf{3}=str;
    StructDate = ComponentStruct(SignalState_cell{Num_zhen});  % 最后一站 最后一帧
    StructDate{1,11}=loc_end;  % 为数据结构体加入定位结果
    targetLocResult = zeros(4,1);
    format1=('InstantPosition-%.3fMHZ-%.3f采样率-%d年%d月%d日-%d时%d分%d秒.mat');
    save('E:\混合数据帧2解析\情景1中间结果存储\Process_M.mat','Process_M');
    save('E:\混合数据帧2解析\情景1中间结果存储\ClusterDoa_cell.mat','ClusterDoa_cell');
    save('E:\混合数据帧2解析\情景1中间结果存储\StateLoc_recoder.mat','StateLoc_recoder');
    save('E:\混合数据帧2解析\情景1中间结果存储\targetLocResult.mat','targetLocResult');
    save(['E:\混合数据帧2解析\CPI时刻的瞬时定位结果\定位情景1\',sprintf(format1,StructDate{2}(1),StructDate{2}(3),StructDate{1}(1),StructDate{1}(2),StructDate{1}(3),StructDate{1}(4),StructDate{1}(5),StructDate{1}(6))],'StructDate')  %%保存该结构体
end
end
   

function [XYZ,flag_isSame]=isSameStation(ByteDate,Num_zhen,Process_M,LBH_ref)
%%%ByteDate 是元胞数组
    LBH_temp=zeros(3,Num_zhen);
    XYZ=[];
    flag_isSame=0;
    if Process_M==0           %%只判断经纬度不换算成坐标
         for i=1:Num_zhen
             date=ByteDate{i};
             jindu_zheng=typecast(date(1,20),'int8');      %%读取经度整数
             jindu_xiaoshu=typecast(date(1,21:22),'int16');%%读取经度小数
             widu_zheng=typecast(date(1,23),'int8');       %%读取纬度整数
             widu_xiaoshu=typecast(date(1,24:25),'int16'); %%读取纬度小数
             Gaodu =typecast(date(1,26:27),'int16');       %%读取高度整数
             jindu=jindu_zheng+(jindu_xiaoshu/2^16);
             widu=widu_zheng+(widu_xiaoshu/2^16);
             LBH_temp(1,i)=jindu;
             LBH_temp(2,i)=widu;
             LBH_temp(3,i)=Gaodu;
         end
             XYZ=LBH_temp(:,2);                       %%此时不做是否同一站判断直接返回坐标 返回中间帧坐标
    else
          for i=1:Num_zhen
             date=ByteDate{i};
             jindu_zheng=typecast(date(1,20),'int8');      %%读取经度整数
             jindu_xiaoshu=typecast(date(1,21:22),'int16');%%读取经度小数
             widu_zheng=typecast(date(1,23),'int8');       %%读取纬度整数
             widu_xiaoshu=typecast(date(1,24:25),'int16'); %%读取纬度小数
             Gaodu =typecast(date(1,26:27),'int16');       %%读取高度整数
             jindu=jindu_zheng+(jindu_xiaoshu/2^16);
             widu=widu_zheng+(widu_xiaoshu/2^16);
             LBH_temp(1,i)=jindu;
             LBH_temp(2,i)=widu;  
             LBH_temp(3,i)=Gaodu;
          end
         
         XYZ_temp=LBM_XYZ(coordinate_LBH,LBH_ref);        %%解算站址坐标
         XYZ_error=XYZ_temp(1:2,2:3)-XYZ_temp(1:2,1);
         error=sqrt(sum(XYZ_error.^2));  count_youxiao=0;
         for i=1:Num_zhen-1
             if error(i)<5.5
                 count_youxiao=count_youxiao+1;
             end
         end
         if count_youxiao==Num_zhen-1
             flag_isSame=1;
             XYZ=XYZ_temp(:,1);                             %%返回站位坐标
         end
    end
end

function CPI_DOAcell=AcquireDoA(SignalState_decoder,zhen_Num)
  % SignalState_decoder 只包含某一个站下的混合数据帧
  CPI_DOAcell=cell(1,zhen_Num);
  for i=1:zhen_Num
     data=SignalState_decoder{i};
     DoA_num=data{41};                        
     CPI_DOA=zeros(DoA_num,3);                
     count_doa=0;
     index=1;
     for j=1:DoA_num
         fangwei=data{41+index};                                % 0-360
          fuyang=data{41+index+1};                              % 0-90
         if fangwei~=0
             CPI_DOA(count_doa+1,1)=fangwei;
             CPI_DOA(count_doa+1,2)=fuyang;     %% 无相对强度数据 暂时填0
             count_doa=count_doa+1;
         end
         index = index+2;
         if (count_doa)>=DoA_num
             break;
         end
     end
     CPI_DOA=CPI_DOA';                         %% 由行数据转换位列数据
     if i==1                                   %% 删除第一帧中的近似角度
         for j=1:DoA_num
            doa1=CPI_DOA(1,j);
            k=j+1;
             while(k<=DoA_num)
                 doa2=CPI_DOA(1,k);
                 if abs(doa1-doa2)<=3
                     CPI_DOA(:,j)=[];
                     DoA_num=DoA_num-1;
                     k=k-1;
                 end
                 k=k+1;
             end
         end
     end
     CPI_DOAcell{1,i}=CPI_DOA;
  end
end

function is=isvaild_n(station_xyz,ref_xyz)
% 情景一未施放升空散射体的瞬时定位点迹计算_第n组有效站位实时判断
% station_xyz各个站位空间坐标，数据类型：当前站位坐标
%ref_xyz 当前第n组测量下的参考站
% 有效判断原则：当前站与上一站相比若距离大于20m 则判断为有效 is=0无效 is=1 有效

    pre=station_xyz-ref_xyz;
    dis=sqrt(pre'*pre);     %
    if dis>=20
        is=1;
    else
        is=0;
    end
end

function DOA_array=DOA_cluster(DOA_cell,I)
% 情景一未施放升空散射体的瞬时定位点迹计算_方位角聚类算法
% 该函数用于UCA阵列在某个站位处DOA聚类，主要目的是过滤非目标的来波方向
% DOA_cell数据类型：元胞矩阵，内容为从I个帧提取出来的DOA数据，每一个cell代表的是一帧的DOA数据第一行方向角，第二行是俯仰角
% I是帧的数量,DOA的单位是度

shape=size(DOA_cell);  %获取元胞数量，正常应该与I相同
DOA_array=zeros(4,100);  %假设获得的角度不会超过100个

cluster_cen=DOA_cell{1};               % 第一个doa矩阵为初始参考聚类中心
cluster_shape=size(cluster_cen);       %获取聚类数量大小
cluster_col=cluster_shape(2);          %记录初始的中心点数量
cluster_cell=cell(1,cluster_col*20);   %创建空的聚类空间

for c=1:cluster_col
    cluster_cell{c}=cluster_cen(:,c);   %依据中心点创建聚类空间 
end
dfai=ones(1,100);                     %作为每次差值的记录,估计中心值不会超过100个
dfai=dfai*500;                        %保证初始元素最大

for c=2:shape(2)
    %从DOA_cell中提取元胞数据用于聚类
    doa_cell=DOA_cell{c};dim=size(doa_cell);                %doa_cell一帧的DOA数据
    for c2=1:dim(2)
        for c3=1:cluster_col
            dfai(c3)=abs(doa_cell(1,c2)-cluster_cen(1,c3));%记录当前doa与中心点的各个差值
        end
          [res,index]=min(dfai);
            if res<=3
                cluster_cell{index}(:,end+1)=doa_cell(:,c2);  %根据最小索引添加到相应聚类空间
                r=size(cluster_cell{index});                %得到index指向的聚类空间
                cluster_cen(:,index)=((r(2)-1)*cluster_cen(:,index)+doa_cell(:,c2))/r(2);
            else
                cluster_cen(:,end+1)=doa_cell(:,c2);  %res是最小值若该值大于3°则作为新的分类中心
                cluster_col=cluster_col+1;            %更新中心点数量
                cluster_cell{cluster_col}=doa_cell(:,c2);   %在聚类空间里面重新创建一个类
            end
            dfai=ones(1,100);                     %重新初始化dfai
            dfai=dfai*500;                       
    end
end

%%%对聚类空间进行检查，过滤掉无效的类
%         cluster_cell(i)=[];    %删除相应的聚类空间
%         cluster_cen(:,i)=[];   %删除相应的中心值   数据正常则中心值和聚类空间值相同
shape2=size(cluster_cell);
col=1;
for i=1:cluster_col
    cell_shape=size(cluster_cell{i});
    if cell_shape(2)>=(I/2)
       cluster_space=cluster_cell{i};
       DOA_array(1,col)=cluster_cen(1,i);        %先放置该方向的中心支
       DOA_array(2,col)=cluster_cen(2,i);        %先放置该方向的中心支
       doa_euler=cluster_space-cluster_cen(:,i); 
       fai_euler=sqrt(sum(doa_euler(1,:).^2));
       theta_euler=sqrt(sum(doa_euler(2,:).^2));
       DOA_array(3,col)=fai_euler;                 %放置该方向的俯仰角均值
       DOA_array(4,col)=theta_euler;               %放置改方向的方向角均值
       col=col+1;
    end
end
DOA_array(:,all(DOA_array==0,1))=[];              %除去全零列
%%%若I较大可以，先进行过滤，从而知道对DOA_CELL预先分配多大空间
end
% 从解析帧中构造结构体
function StructDate = ComponentStruct(Sigframe_decoded)
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
 P_num = Sigframe_decoded{41};
 StructDate{8} = zeros(8,3);     % DOA信息 8*3 第一行 DOA的数量
                                % 方向角 俯仰角 来波方向强度
StructDate{8}(1,1)=P_num;
StructDate{9}=zeros(5,5);                   % zeros 5*5  T散射体信息
refdoa_order = Sigframe_decoded{65};
StructDate{10}=zeros(1,3);                  % 参考通道信号来波方向 DOA信息中的字段
StructDate{10}(1:2) = StructDate{8}(refdoa_order+1,1:2);
StructDate{10}(3) =  refdoa_order;%
end

function location=line_location(station_n,DOA_n)
% 情景一未施放升空散射体的瞬时定位点迹计算_双站二维射线交汇定位
% station_n 第n组有效的双站坐标x、y、z
% DOA_n数据类型：元胞矩阵，内容为第n组双站各自的DOA数据 cell：方向角、俯仰角、方向角均方误差、俯仰角均方误差
% location 应该包括交汇定位结果，误差信任度，俯仰角
% 该函数用于第n组前后两站直线交汇定位，并记录记录定位点坐标和误差信任度
% 注意： 送入的角度必须保证0-360

%%%初始化
state_1=station_n(1:2,1);  %%创建两个站的独立坐标
state_2=station_n(1:2,2);
DOA_1=DOA_n{1};          %%存储两个站的DOA数据
DOA_2=DOA_n{2};
shape1=size(DOA_1);      %%记录大小
shape2=size(DOA_2);
loc_cell=zeros(1,5);     %%交汇定位的最小单元 （方向角1，方向角2,均值方差1，均值方差2;俯仰角）
location=zeros(4,min(shape1(2),shape2(2)));  Index=1;
for col=1:shape1(2)
    loc_cell(1)=DOA_1(1,col);  %提取一个最小单元
    loc_cell(3)=DOA_1(3,col);
    for col2=1:shape2(2)   
        loc_cell(2)=DOA_2(1,col2);
        loc_cell(4)=DOA_2(3,col2);
        loc_cell(5)=DOA_2(2,col2);  %提取俯仰角
        if (abs(loc_cell(1)-loc_cell(2))<=1)   %%方案中判断条件是相等，这里改为小于1度即认作相同
            continue;
        end
        Km_1=tan(loc_cell(1)*(pi/180));
          Km=tan(loc_cell(2)*(pi/180));
        x_j=(state_2(2)-state_1(2)-Km*state_2(1)+Km_1*state_1(1))/(Km_1-Km);
        y_j=(Km_1*state_2(2)-Km*state_1(2)-Km*Km_1*(state_2(1)-state_1(1)))/(Km_1-Km);
        if (loc_cell(2)>=0 && loc_cell(2)<90)
            if (x_j>state_2(1) && y_j>=state_2(2))
                flag=1 ;                     %%标志位，用于判断是否信任定位点 1为信任
            else
                flag=0;
            end
        elseif (loc_cell(2)>=90 && loc_cell(2)<180)
            if (x_j<=state_2(1) && y_j>state_2(2))
                flag=1 ;                     %%标志位，用于判断是否信任定位点 1为信任
            else
                flag=0;
            end
         elseif (loc_cell(2)>=180 && loc_cell(2)<270)
            if (x_j<state_2(1) && y_j<=state_2(2))
                flag=1 ;                     %%标志位，用于判断是否信任定位点 1为信任
            else
                flag=0;
            end
        elseif (loc_cell(2)>=180 && loc_cell(2)<270)
             if (x_j<state_2(1) && y_j<=state_2(2))
                flag=1 ;                     %%标志位，用于判断是否信任定位点 1为信任
            else
                flag=0;
             end
        elseif (loc_cell(2)>=270 && loc_cell(2)<360)
            if (x_j>=state_2(1) && y_j<state_2(2))
                flag=1 ;                     %%标志位，用于判断是否信任定位点 1为信任
            else
                flag=0;
            end
        end
        if (flag==1)
            %%计算误差信任度函数
            m_1=sec(loc_cell(1)*(pi/180));
            m=sec(loc_cell(2)*(pi/180));
            A_j=(state_2(2)-state_1(2)-(state_2(1)-state_1(1))*Km_1)/(Km_1-Km)^2*m^2;
            B_j=(state_2(2)-state_1(2)-(state_2(1)-state_1(1))*Km)/(Km_1-Km)^2*m_1^2;
            dx=A_j*loc_cell(4)+B_j*loc_cell(3);
            dy=Km_1*A_j*loc_cell(4)+Km*B_j*loc_cell(3);
            S_j=1/(dx*dy);          %%计算误差信任度
            location(1,Index)=x_j;  %%记录有效定位点
            location(2,Index)=y_j;
            location(3,Index)=S_j;
            location(4,Index)=loc_cell(5);
            Index=Index+1;
        end
    end
end
end

function [sta_fw,sta_fy,alph,beta] = AntennaCor2StandardStationCor(H,fw,fy,ant_fw,ant_fy)

% 根据当前天线阵所在的xyz磁场分量、磁偏角、磁倾角信息，将天线阵doa估计得到的角度
% 信息转换为标准站心坐标系下的方位角和俯仰角。后面的定位都是在标准站心坐标系下进行的

% H：三磁场分量，(Hx,Hy,Hz)1*3的分量
% fw：磁偏角，角度，与y轴夹角（地磁N级）
% fy：磁倾角，角度，与xoy面的夹角
% ant_fw：天线阵坐标系中的方位角，角度
% ant_fy：天线阵坐标系中的俯仰角，角度
% 
% sta_fw：站心坐标系的方位角，角度
% sta_fy：站心坐标系的俯仰角，角度

H_value = sqrt(H*H');
x1 = [H_value*cosd(fw)*cosd(fy);H_value*sind(fw)*cosd(fy);H_value*sind(fy)];
x2 = H';
A = (x1*x2')*(x2*x2')^(-1);
beta = acosd(H_value*sind(fy)/sqrt(H(1:2)*H(1:2)')) + atand(H(1)/H(3));
alph = acosd(H_value*cosd(fy)*sind(fw)/sqrt((H(1,1)*cosd(beta) - H(1,3)*sind(beta))^2 + H(1,2)^2)) - atand(H(1,3)/(H(1,1)*cosd(beta) - H(1,3)*sind(beta)));
ant = [cosd(ant_fy)*cosd(ant_fw);cosd(ant_fy)*sind(ant_fw);sind(ant_fy)];
sta = A*ant;
sta_fw = atand(sta(2,1)/sta(1,1));
sta_fy = asind(sta(3,1));


end


function [XYZ,Target_V, trace_cell,vxyz_cell]=Qin_LocByTraceOfPointFormal(Loc_cell,CPI_num,Vmin,Vmax,Vz_min,Vz_max,detax,detay,detaz,detaT,GRADE_MAX,K)
%%%跨CPI时刻点迹分选和滤波平滑
%%%参数说明： Loc_cell: K帧定位结果 元胞矩阵 列是数据  未决定是否存在无效数据
%%%          Vmin：  水平运动最小速度
%%%          Vmax：  水平运动最大速度
%%%          Vz_min：升降运动最小速度
%%%          Vz_max: 升降运动最大速度
%%%           detax: x坐标波门
%%%           detay：y坐标波门
%%%           detaz: z坐标波门
%%%           GRADE_MAX: 目标的初始分数
%%%           K：递推最小二乘点数
%%%          zhen_size: 帧的数量
Time_nums = CPI_num;                              % 记录当前正在处理的CPI时刻，初始送入CPI_num应该等于0
Targets_nums=load('E:\混合数据帧2解析\点迹分选中间结果存储\Targets_nums.mat').Targets_nums;   % 记录当前点迹初始点数量
All_targets =load('E:\混合数据帧2解析\点迹分选中间结果存储\All_targets.mat').All_targets;    % All_targets存放着目标轨迹，轨迹点数、中心波门、分数、轨迹创建CPI_时刻
  zhen_size = K;
for i=1:zhen_size
  temp_loc=Loc_cell{i};
  temp_loc_size=size(temp_loc);
  Point_cpi=zeros(5,temp_loc_size(2));       %% 当前cpi时刻点迹矩阵 一个编号、关联标志flag、坐标值x,y,z
  Point_cpi(3:5,1:end)=temp_loc(1:3,1:end);  %% 建立当前点迹组合 ALL_points
  Point_cpi(1,1:end)=1:temp_loc_size(2);
  if CPI_num ==0,Time_nums=Time_nums+1; end 
  if Time_nums==1
    for j=1:temp_loc_size(2)
       All_targets{j}=cell(1,7);             %%一个编号name、建立时间beginning 计分器grade 坐标计数器nums 点迹编号数组 点迹坐标数组 当前时刻波门中心点坐标
       All_targets{j}{1}=Point_cpi(1,j);     % 编号name
       All_targets{j}{2}=Time_nums;          % 记录当前点迹初始点的建立时间beginning 即从属于第几帧
       All_targets{j}{3}=GRADE_MAX;          % 目标在无点迹匹配存活的最大时间
       All_targets{j}{4}=1;                  % 坐标计数器。即以当前点迹为初始点的轨迹包含的点迹数量
       All_targets{j}{5}=zeros(1,200);       % 编号数组大小
       All_targets{j}{6}=zeros(3,200);       % 记录当前轨迹包含的点迹坐标
       All_targets{j}{6}(:,1)=Point_cpi(3:5,j);
       All_targets{j}{7}=Point_cpi(3:5,j);   % 当前轨迹加入条件：波门中心点，待加入点会计算与波门中心点的dxy
       Targets_nums=Targets_nums+1;
    end
    continue;
  end
   Dxy_mp=ones(Targets_nums,temp_loc_size(2))*9999; % 目标点迹距离矩阵 行数：目标点 列数：点迹点
   for j=1:Targets_nums
       temp_tarcell=All_targets{j};     
       trace_nums=temp_tarcell{4};      % 取出当前轨迹包含的点迹数
       grade=temp_tarcell{3};           % 取出的当前轨迹的分数
       center_point=temp_tarcell{7};    
     for k=1:temp_loc_size(2) 
       point_xyz=temp_loc(:,k);         %%当前cpi第k列点迹坐标 
       dxy=sqrt((center_point(1:2)-point_xyz(1:2))'*(center_point(1:2)-point_xyz(1:2)));
       dz=sqrt((center_point(3)-point_xyz(3))'*(center_point(3)-point_xyz(3)));
       if trace_nums==1
           N=GRADE_MAX;
           if (dxy>Vmin*detaT*(N-grade+1) && dxy<=Vmax*detaT*(N-grade+1)) && (dz>=Vz_min*detaT*(N-grade+1) && dz<=Vz_max*detaT*(N-grade+1))
               Dxy_mp(j,k)=dxy;
           end
       else
           N=min(trace_nums,K);
           x=sqrt((center_point(1)-point_xyz(1))'*(center_point(1)-point_xyz(1)));
           y=sqrt((center_point(2)-point_xyz(2))'*(center_point(2)-point_xyz(2)));
           z=sqrt((center_point(3)-point_xyz(3))'*(center_point(3)-point_xyz(3)));
           if (x<=detax*N/grade) && (y<=detay*N/grade) && (z<=detaz*N/grade)
               Dxy_mp(j,k)=dxy;
           end
       end
     end
   end
   t=zeros(1,Targets_nums);   %%索引是目标 内容是点迹
   p=zeros(1,size(Dxy_mp,2));   %%索引是点迹 内容是目标
   for j=1:Targets_nums       %%用t记录点迹
       [value,index]=min(Dxy_mp(j,:));
       if value~=9999
       t(1,j)=index;
       end
     if (t(1,j)~=0)
       if p(1,t(1,j))==0      %%为p存放合适目标点
           p(1,t(1,j))=j;
       else
           dxy1=Dxy_mp(j,t(1,j));
           dxy2=Dxy_mp(p(1,t(1,j)),t(1,j));
           grade1=All_targets{j}{3};
           grade2=All_targets{p(1,t(1,j))}{3};
           if dxy1<dxy2    %%竞争成功
             Dxy_mp(p(1,t(1,j)),t(1,j))=9999;
             p(1,t(1,j))=j;
           elseif dxy1==dxy2
               if grade1 >grade2
                   Dxy_mp(p(1,t(1,j)),t(1,j))=9999;
                   p(1,t(1,j))=j;
               elseif grade1==grade2
                    name1=All_targets{j}{1};
                    name2=All_targets{p(1,t(1,j))}{1};
                   if name1<name2
                     Dxy_mp(p(1,t(1,j)),t(1,j))=9999;
                     p(1,t(1,j))=j;
                   end
               end
           else       %%竞争失败
               Dxy_mp(j,t(1,j))=9999;
           end
       end
     end
   end
   t2=zeros(2,Targets_nums);
   for j=1:Targets_nums
       if t(1,j)~=0
           if p(1,t(1,j))==j
               t2(1,j)=t(1,j);
               t2(2,j)=1;
           else
               [dxy,index]=min(Dxy_mp(j,:)) ;  % 当前目标点的点迹已经被取代，从Dxy寻找下一个
               if dxy~=9999
                t2(1,j)=index;
                t2(2,j)=1;   
               end
           end
       end
   end
   % 开始更新All_targets
   j=1; 
   while j<=Targets_nums     
      if t2(1,j)~=0
         Point_cpi(2,t2(1,j))=1;
         trace_nums=All_targets{j}{4};
         beginning=All_targets{j}{2};
         All_targets{j}{5}(1,Time_nums-beginning+1)=t2(1,j);
         All_targets{j}{6}(:,Time_nums-beginning+1)= Point_cpi(3:5,t2(1,j));
         if (Time_nums-beginning+1-trace_nums)>=2
             trace_nums=trace_nums+1;
             All_targets{j}{4}=All_targets{j}{4}+1;
             All_targets{j}{6}(:,trace_nums)= (All_targets{j}{6}(:,1)+All_targets{j}{6}(:,Time_nums-beginning+1))*0.5;
         end
         All_targets{j}{4}=All_targets{j}{4}+1;
         trace_nums=All_targets{j}{4};
         N=min(trace_nums,K);
         center_point=waituichazhi(All_targets{j}{6}(:,trace_nums-N+1:trace_nums),detaT,N); %%利用外推差值计算波门中心点

         All_targets{j}{7}=center_point;
         All_targets{j}{3}=All_targets{j}{3}+1;
         if All_targets{j}{3}>GRADE_MAX
             All_targets{j}{3}=GRADE_MAX;
         end
      elseif t2(1,j)==0
         trace_nums=All_targets{j}{4};
         if trace_nums==1
            All_targets{j}{3}=All_targets{j}{3}-1; 
          if All_targets{j}{3}==0
              All_targets(j)=[];     %%删除第j个目标
              t2(:,j)=[];
              j=j-1;
              Targets_nums=Targets_nums-1;
          end
         elseif trace_nums>1               %%nums不为1 则用波门点代替下一个点迹
            All_targets{j}{6}(:,trace_nums+1)= All_targets{j}{7};
            N=min(trace_nums+1,K);
            trace_nums = trace_nums+1;
            center_point=waituichazhi(All_targets{j}{6}(:,trace_nums-N+1:trace_nums),detaT,N); %%利用外推差值计算波门中心点
            All_targets{j}{7}=center_point;
            All_targets{j}{4}=All_targets{j}{4}+1;
            All_targets{j}{3}=All_targets{j}{3}-1; 
          if All_targets{j}{3}==0
              All_targets(j)=[];     %%删除第j个航迹
              t2(:,j)=[];
              j=j-1;
              Targets_nums=Targets_nums-1;
          end
         end
      end
      j=j+1;
   end
   flag=Point_cpi(2,:);                                   %%将孤立点重新作为目标点
   for j=1:temp_loc_size(2)
       if flag(j)==0
       All_targets{Targets_nums+1}=cell(1,7);             %%一个编号name、建立时间beginning 计分器grade 坐标计数器nums 点迹编号数组 点迹坐标数组 当前时刻波门中心点坐标
       All_targets{Targets_nums+1}{1}=Point_cpi(1,j);     %%编号name
       All_targets{Targets_nums+1}{2}=Time_nums;          %%一个建立时间beginning
       All_targets{Targets_nums+1}{3}=GRADE_MAX;          %%一个计分器  目标在无点迹匹配存活的最大时间
       All_targets{Targets_nums+1}{4}=1;                  %%坐标计数器
       All_targets{Targets_nums+1}{5}=zeros(1,200);       %%编号数组大小
       All_targets{Targets_nums+1}{6}=zeros(3,200);       %%坐标数组大小
       All_targets{Targets_nums+1}{6}(:,1)=Point_cpi(3:5,j);
       All_targets{Targets_nums+1}{7}=Point_cpi(3:5,j);   %%当前时刻波门中心点坐标
       Targets_nums=Targets_nums+1;   
       end
   end
   if Time_nums>=zhen_size     %%读取的帧数大于等于K值则
       break;
   end
end
save('E:\混合数据帧2解析\点迹分选中间结果存储\All_targets.mat','All_targets');
save('E:\混合数据帧2解析\点迹分选中间结果存储\Targets_nums.mat','Targets_nums');
% 利用反向递推进行最终定位点
XYZ_cell=cell(1,Targets_nums);
Target_V_cell=cell(1,Targets_nums);
count=0;
    for i=1:Targets_nums
        targets=All_targets{i};
           trace_nums=targets{4};
        if trace_nums>=(zhen_size-1)
            count=count+1;
            trace_point=targets{6};
            [xyz,V_xyz]=fanxiangditui(trace_point,trace_nums,K,detaT,Time_nums);
            XYZ_cell{1,count}=xyz;
            Target_V_cell{1,count}=V_xyz;    
        end    
    end
  XYZ = zeros(3,200);
  Target_V=zeros(3,200);
  index=1; count=1;
  trace_cell={};
  vxyz_cell = {};
  for i=1:Targets_nums
      xyz = XYZ_cell{i};
      V_xyz = Target_V_cell{1,i};
      columns = size(xyz,2);
      if size(xyz,2)>0
          XYZ(:,index:index+columns-1) = xyz;
          Target_V(:,index:index+columns-1)=V_xyz;
          index = index+columns;
          trace_cell{count}= XYZ_cell{i};
          vxyz_cell{count} = Target_V_cell{i};
          count =count+1;
      end
  end
  XYZ(:,all(XYZ==0,1))=[];
  Target_V(:,all(Target_V==0,1))=[];
end

function center_point=waituichazhi(loc,detaT,N)
     center_point=zeros(3,1);
     x=loc(1,:)';
     y=loc(2,:)'; 
     z=loc(3,:)';
     H=ones(N,2);  H(:,2)=-flipud((1:N)')*detaT;
     x2=(H'*H)\H'*x;
     y2=(H'*H)\H'*y;
     z2=(H'*H)\H'*z;
     center_point(1,1)=x2(1,1);
     center_point(2,1)=y2(1,1);
     center_point(3,1)=z2(1,1);
end

function [xyz,V_xyz]=fanxiangditui(trace_point,nums,K,detaT,Time_nums)
     H = zeros(K,2);
     H(1:K,1) = ones(K,1);
     H(2:K,2) = detaT*(1:K-1)';
     xyz = zeros(3,nums);
     V_xyz = zeros(3,nums);
     i=Time_nums;
     x=trace_point(1,i-K+1:i)';
     y=trace_point(2,i-K+1:i)';
     z=trace_point(3,i-K+1:i)';
     x2=(H'*H)\H'*x;
     y2=(H'*H)\H'*y;
     z2=(H'*H)\H'*z;
     xyz(1,i-K+1)=x2(1,1);    V_xyz(1,i-K+1)=x2(2,1);
     xyz(2,i-K+1)=y2(1,1);    V_xyz(2,i-K+1)=y2(2,1);
     xyz(3,i-K+1)=z2(1,1);    V_xyz(3,i-K+1)=z2(2,1);
end

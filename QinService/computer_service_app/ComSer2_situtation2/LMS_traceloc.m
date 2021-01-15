function [loc_shunshi,couple_array]=LMS_traceloc(point_cluster,S_xyz,S_varray,f_array,doa_ref,flag_caiyang,flag_yundong)
 %%%应用于情景2 释放升空散射体但未检测到散射体，最小二乘法进行TDOA-DOA定位  FDOA-DOA定位
 %%%point_cluster 是Q个散射体得到的所有可能定位点迹,胞元矩阵 ，行是数据
 %%%行是组合结果：TDOA时差、等效信噪比_tdoa、FDOA时差、等效信噪比_fdoa、方向角_ref，俯仰角_ref 散射体索引
 %%%S_xyz 是转换坐标后的散射体坐标， 有效+无效  列是数据
 %%%S_varray 散射体的速度， 列是数据
 %%%doa_ref：1*2 参考通道doa数据   fla_caiyang:  1: 高采样   0： 低采样   
 %%%flag_yundong 0: 不计算速度 1：计算速度
 %%%loc_shunshi 定位矩阵点， couple_array 有效的定位组合 元胞矩阵 列是数据
 
 shape=size(point_cluster);
 num=shape(2);                %%%得到胞元的数量

 if num==1
     Index=-1;
      if flag_caiyang==1          %%%分为高采样 低采样 不同定位形式
         si_box=point_cluster{1}; 
       if si_box(1,1)~=-1              %%%检查散射体在第二次筛选后的有效性
         tao_array=si_box(:,1);
         tao_array(:,2)=si_box(:,7);
         s_xyz=S_xyz;   %%%提取出相应散射体坐标
         [loc_shun,Index]=S1_loc_Gaocai(s_xyz,tao_array,doa_ref);   %%%返回瞬时定位， 和配对索引
       else
         loc_shun=-1; 
       end 
      else
         si_box=point_cluster{1};
         if size(si_box,1)~=0
             fd_array=si_box(:,3);
             fd_array(:,2)=si_box(:,7);
             s_xyz=S_xyz;   %%%将所有散射体送入
             f=f_array;
             s_v=S_varray;
             [loc_shun,Index]=S1_loc_Dicaij(s_xyz,s_v,fd_array,doa_ref,f);   %%%返回瞬时定位， 和配对索引
         else
          loc_shun=-1;     %%%%区别于=0的情况  -1 代表LSM没有找到有效定位点
         end
      end 

     if loc_shun(1,1)~=-1
         value2=size(Index);
         loc_shunshi=loc_shun;
         couple_array=cell(1,1);     %%%组合数据收集   列是数据
         couple_array{1}=zeros(7,value2(2));
         for i=1:value2(2)
            couple_array{1}(:,i)=si_box(Index(i),:)';   %%% 记录求解出的定位点所利用的组合参数
         end
     else
     loc_shunshi=-1;
     couple_array=-1;
     end
 end
 
 if num==2
     si_box=point_cluster{1};
     sisi_box=point_cluster{2};
     loc_shun1=-1;  loc_shun2=-1;
   if si_box(1,1) ~=-1
     if flag_caiyang==1       %%%分为高采样 低采样 不同定位形式
         tao_array=si_box(:,1);
         tao_array(:,2)=si_box(:,7);
         s_xyz=S_xyz;   %%%将所有散射体送入
         [loc_shun1,Index]=S1_loc_Gaocai(s_xyz,tao_array,doa_ref);   %%%返回瞬时定位， 和配对索引   
     else
         fd_array=si_box(:,3);
         fd_array(:,2)=si_box(:,7);
         s_xyz=S_xyz;   %%%将所有散射体送入
         f=f_array;
         s_v=S_varray;
         if (flag_yundong==0)
         [loc_shun1,Index]=S1_loc_Dicaij(s_xyz,s_v,fd_array,doa_ref,f);   %%%返回瞬时定位， 和配对索引
         else
          %%%%%运动目标散射体
         end
     end
   else
      loc_shun1=-1;
   end
   
   if sisi_box(1,1)~=-1
       if flag_caiyang==1       %%%分为高采样 低采样 不同定位形式
         tao_array=sisi_box(:,1);
         tao_array(:,2)=sisi_box(:,5);
         tao_array(:,3:4)=sisi_box(:,11:12);
         s_xyz=S_xyz;   %%%将所有散射体送入
         [loc_shun2,Index2]=S1S2_loc_Gaocai(s_xyz,tao_array,doa_ref);   %%%返回瞬时定位， 和配对索引   
      else
         fd_array=sisi_box(:,3);
         fd_array(:,2)=sisi_box(:,7);
         fd_array(:,3:4)=sisi_box(:,11:12);
         s_xyz=S_xyz;   %%%将所有散射体送入
         f=f_array;
         s_v=S_varray;
         [loc_shun2,Index2]=S1S2_loc_Dicaij(s_xyz,s_v,fd_array,doa_ref,f);   %%%返回瞬时定位， 和配对索引
      end
   else
       loc_shun2=-1;
   end
   
    couple_array=cell(1,2);              %%组合数据收集
    couple_array{1}=0;  
    couple_array{2}=0; 
   if loc_shun1(1,1)~=-1 && loc_shun2(1,1)~=-1
     value1=size(Index);   value2=size(Index2);
     loc_shunshi=[loc_shun1 loc_shun2];
     couple_array{1}=zeros(7,value1(2));
     couple_array{2}=zeros(12,value2(2));
         for i=1:value1(2)
            couple_array{1}(:,i)=si_box(Index(i),:)';
         end
         for i=1:value2(2)
            couple_array{2}(:,i)=sisi_box(Index2(i),:)';
         end
   elseif loc_shun1(1,1)~=-1
     value1=size(Index);  
       couple_array{1}=zeros(7,value1(2));
       value1=size(Index); 
       loc_shunshi=loc_shun1;
         for i=1:value1(2)
            couple_array{1}(:,i)=si_box(Index(i),:)';
         end
   elseif loc_shun2(1,1)~=-1
       value2=size(Index2);
       loc_shunshi=loc_shun2;
       couple_array{2}=zeros(12,value2(2));
         for i=1:value2(2)
            couple_array{2}(:,i)=sisi_box(Index2(i),:)';
         end
   else
       loc_shunshi=-1;
       couple_array=-1;
   end
     
 end
 
 if num==3
     si_box=point_cluster{1};                   %% 行数是数据组合数 列是数据
     sisi_box=point_cluster{2};
     sisisi_box=point_cluster{3};
 
     
   if si_box(1,1) ~=-1
     if flag_caiyang==1       %%%分为高采样 低采样 不同定位形式
         tao_array=si_box(:,1);
         tao_array(:,2)=si_box(:,7);
         s_xyz=S_xyz;   %%%将所有散射体送入
         [loc_shun1,Index1]=S1_loc_Gaocai(s_xyz,tao_array,doa_ref);   %%%返回瞬时定位， 和配对索引   
     else
         fd_array=si_box(:,3);
         fd_array(:,2)=si_box(:,7);
         s_xyz=S_xyz;   %%%将所有散射体送入
         f=f_array;
         s_v=S_varray;
         [loc_shun1,Index1]=S1_loc_Dicaij(s_xyz,s_v,fd_array,doa_ref,f);   %%%返回瞬时定位， 和配对索引
     end
   else
      loc_shun1=-1;
   end
   
   if sisi_box(1,1)~=-1
       if flag_caiyang==1       %%%分为高采样 低采样 不同定位形式
         tao_array=sisi_box(:,1);
         tao_array(:,2)=sisi_box(:,5);
         tao_array(:,3:4)=sisi_box(:,11:12);
         s_xyz=S_xyz;   %%%将所有散射体送入
         [loc_shun2,Index2]=S1S2_loc_Gaocai(s_xyz,tao_array,doa_ref);   %%%返回瞬时定位， 和配对索引   
      else
         fd_array=sisi_box(:,3);
         fd_array(:,2)=sisi_box(:,7);
         fd_array(:,3:4)=sisi_box(:,11:12);
         s_xyz=S_xyz;   %%%将所有散射体送入
         f=f_array;
         s_v=S_varray;
         [loc_shun2,Index2]=S1S2_loc_Dicaij(s_xyz,s_v,fd_array,doa_ref,f);   %%%返回瞬时定位， 和配对索引
      end
   else
       loc_shun2=-1;
   end

    if sisisi_box(1,1)~=-1
       if flag_caiyang==1       %%%分为高采样 低采样 不同定位形式
         tao_array=sisisi_box(:,1);
         tao_array(:,2)=sisisi_box(:,5);
         tao_array(:,3)=sisisi_box(:,9);
         tao_array(:,4:6)=sisisi_box(:,15:17);
         s_xyz=S_xyz;   %%%将所有散射体送入
         [loc_shun3,Index3]=S1S2S3_loc_Gaocai(s_xyz,tao_array,doa_ref);   %%%返回瞬时定位， 和配对索引   
      else
         fd_array=sisisi_box(:,3);
         fd_array(:,2)=sisisi_box(:,7);
         fd_array(:,3)=sisisi_box(:,11);
         fd_array(:,4:6)=sisisi_box(:,15:17);
         s_xyz=S_xyz;   %%%将所有散射体送入
         f=f_array;
         s_v=S_varray;
         [loc_shun3,Index3]=S1S2S3_loc_Dicaij(s_xyz,s_v,fd_array,doa_ref,f);   %%%返回瞬时定位， 和配对索引
      end
   else
       loc_shun3=-1;
   end

loc_shunshi=zeros(3,1);   
couple_array=cell(1,3);
couple_array{1}=0;   
couple_array{2}=0;  
couple_array{3}=0;
if loc_shun1(1,1)~=-1
    value1=size(Index1); 
   loc_shunshi=[loc_shunshi loc_shun1];
   couple_array{1}=zeros(7,value1(2));     %%收集数据
     for i=1:value1(2)
        couple_array{1}(:,i)=si_box(Index1(i),:)';
     end
end
if loc_shun2(1,1)~=-1
   value2=size(Index2); 
   loc_shunshi=[loc_shunshi loc_shun2];
      couple_array{2}=zeros(12,value2(2));     %%收集数据
     for i=1:value2(2)
        couple_array{2}(:,i)=sisi_box(Index2(i),:)';
     end
end
if loc_shun3(1,1)~=-1
    value3=size(Index3); 
   loc_shunshi=[loc_shunshi loc_shun3];
   couple_array{3}=zeros(17,value3(2));     %%收集数据
     for i=1:value3(2)
        couple_array{3}(:,i)=sisisi_box(Index3(i),:)';
     end
end
col=size(loc_shunshi);
if col(2)>1
  loc_shunshi(:,1)=[];
else
  loc_shunshi=-1;
  couple_array=-1;
end

 end
 
 if num==4
     si_box=point_cluster{1};
     sisi_box=point_cluster{2};
     sisisi_box=point_cluster{3};
     sisisisi_box=point_cluster{4};
    
   if si_box(1,1) ~=-1
     if flag_caiyang==1       %%%分为高采样 低采样 不同定位形式
         tao_array=si_box(:,1);
         tao_array(:,2)=si_box(:,7);
         s_xyz=S_xyz;   %%%将所有散射体送入
         [loc_shun1,Index1]=S1_loc_Gaocai(s_xyz,tao_array,doa_ref);   %%%返回瞬时定位， 和配对索引   
     else
         fd_array=si_box(:,3);
         fd_array(:,2)=si_box(:,7);
         s_xyz=S_xyz;   %%%将所有散射体送入
         f=f_array;
         s_v=S_varray;
         [loc_shun1,Index1]=S1_loc_Dicaij(s_xyz,s_v,fd_array,doa_ref,f);   %%%返回瞬时定位， 和配对索引
     end
   else
      loc_shun1=-1;
   end
   
   if sisi_box(1,1)~=-1
       if flag_caiyang==1       %%%分为高采样 低采样 不同定位形式
         tao_array=sisi_box(:,1);
         tao_array(:,2)=sisi_box(:,5);
         tao_array(:,3:4)=sisi_box(:,11:12);
         s_xyz=S_xyz;   %%%将所有散射体送入
         [loc_shun2,Index2]=S1S2_loc_Gaocai(s_xyz,tao_array,doa_ref);   %%%返回瞬时定位， 和配对索引   
      else
         fd_array=sisi_box(:,3);
         fd_array(:,2)=sisi_box(:,7);
         fd_array(:,3:4)=sisi_box(:,11:12);
         s_xyz=S_xyz;   %%%将所有散射体送入
         f=f_array;
         s_v=S_varray;
         [loc_shun2,Index2]=S1S2_loc_Dicaij(s_xyz,s_v,fd_array,doa_ref,f);   %%%返回瞬时定位， 和配对索引
      end
   else
       loc_shun2=-1;
   end

    if sisisi_box(1,1)~=-1
       if flag_caiyang==1       %%%分为高采样 低采样 不同定位形式
         tao_array=sisisi_box(:,1);
         tao_array(:,2)=sisisi_box(:,5);
         tao_array(:,3)=sisisi_box(:,9);
         tao_array(:,4:6)=sisisi_box(:,15:17);
         s_xyz=S_xyz;   %%%将所有散射体送入
         [loc_shun3,Index3]=S1S2S3_loc_Gaocai(s_xyz,tao_array,doa_ref);   %%%返回瞬时定位， 和配对索引   
      else
         fd_array=sisisi_box(:,3);
         fd_array(:,2)=sisisi_box(:,7);
         fd_array(:,3)=sisisi_box(:,11);
         fd_array(:,4:6)=sisisi_box(:,15:17);
         s_xyz=S_xyz;   %%%将所有散射体送入
         f=f_array;
         s_v=S_varray;
         [loc_shun3,Index3]=S1S2S3_loc_Dicaij(s_xyz,s_v,fd_array,doa_ref,f);   %%%返回瞬时定位， 和配对索引
      end
   else
       loc_shun3=-1;
   end

  if sisisisi_box(1,1)~=-1
       if flag_caiyang==1       %%%分为高采样 低采样 不同定位形式
         tao_array=sisisi_box(:,1);
         tao_array(:,2)=sisisi_box(:,5);
         tao_array(:,3)=sisisi_box(:,9);
         tao_array(:,4)=sisisi_box(:,13);
         s_xyz=S_xyz;   %%%将所有散射体送入
         [loc_shun4,Index4]=S1S2S3S4_loc_Gaocai(s_xyz,tao_array,doa_ref);   %%%返回瞬时定位， 和配对索引   
      else
         fd_array=sisisi_box(:,3);
         fd_array(:,2)=sisisi_box(:,7);
         fd_array(:,3)=sisisi_box(:,11);
         fd_array(:,4)=sisisi_box(:,15);
         s_xyz=S_xyz;   %%%将所有散射体送入
         f=f_array;
         s_v=S_varray;
         [loc_shun4,Index4]=S1S2S3S4_loc_Dicaij(s_xyz,s_v,fd_array,doa_ref,f);   %%%返回瞬时定位， 和配对索引
      end
   else
       loc_shun4=-1;
   end

loc_shunshi=zeros(3,1);
couple_array=cell(1,4);
couple_array{1}=0;   couple_array{2}=0;  couple_array{3}=0; couple_array{4}=0;
if loc_shun1(1,1)~=-1
    value1=size(Index1); 
   loc_shunshi=[loc_shunshi loc_shun1];
      couple_array{1}=zeros(7,value1(2));     %%收集数据
     for i=1:value1(2)
        couple_array{1}(:,i)=si_box(Index1(i),:);
     end
end
if loc_shun2(1,1)~=-1
   value2=size(Index2); 
   loc_shunshi=[loc_shunshi loc_shun2];
      couple_array{2}=zeros(12,value2(2));     %%收集数据
     for i=1:value2(2)
        couple_array{2}(:,i)=sisi_box(Index2(i),:)';
     end
end
if loc_shun3(1,1)~=-1
   value3=size(Index3); 
   loc_shunshi=[loc_shunshi loc_shun3];
   couple_array{3}=zeros(17,value3(2));     %%收集数据
     for i=1:value3(2)
        couple_array{3}(:,i)=sisisi_box(Index3(i),:)';
     end
end
if loc_shun4(1,1)~=-1
   value4=size(Index4); 
   loc_shunshi=[loc_shunshi loc_shun4];
      couple_array{4}=zeros(18,value4(2));     %%收集数据
     for i=1:value4(2)
        couple_array{4}(:,i)=sisisisi_box(Index4(i),1:18)';
     end
end
col=size(loc_shunshi);
if col(2)>1
  loc_shunshi(:,1)=[];
else
  loc_shunshi=-1;
  couple_array=-1;
end
     
end
 
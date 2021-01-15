function [loc_shunshi,couple_array]=Qin3_LMStraceloc(point_cluster,S_xyz,S_ref,S_varray,S_vref,f_array,doa_ref,flag_caiyang,flag_branch,Q,T)
 %%%应用于情景2 释放升空散射体但未检测到散射体，最小二乘法进行TDOA-DOA定位  FDOA-DOA定位
 %%%point_cluster 是Q个散射体得到的所有可能定位点迹,胞元矩阵 ，行是数据
 %%%行是组合结果：TDOA时差、等效信噪比_tdoa、FDOA时差、等效信噪比_fdoa、方向角_ref，俯仰角_ref 散射体索引
 %%%S_xyz 是转换坐标后的散射体坐标， 有效+无效  列是数据
 %%%S_varray 散射体的速度， 列是数据
 %%%doa_ref：1*2 参考通道doa数据   fla_caiyang:  1: 高采样   0： 低采样   
 %%%flag_yundong 0: 不计算速度 1：计算速度
 %%%loc_shunshi 定位矩阵点， couple_array 有效的定位组合 元胞矩阵 列是数据
 
 %%%%%%%%%%分支1%%%%%%%%%%%%
 if flag_branch==1
     branch1_box=point_cluster{1};       %%%在分支一胞元矩阵只包含一个矩阵
    if branch1_box(1,1) ~=-1
       doa_rd=branch1_box(:,5:6);
     if flag_caiyang==1                  %%%分为高采样 低采样 不同定位形式
         tao_array=branch1_box(:,1);
         tao_array(:,2)=branch1_box(:,7);
         [loc_shun1,Index1]=Qin3_branch1_locGaocai(S_ref(:,1),tao_array,doa_rd);   %%%返回瞬时定位， 和配对索引   
     else
         fd_array=branch1_box(:,3);
         fd_array(:,2)=branch1_box(:,7);
         f=f_array;
         [loc_shun1,Index1]=Qin3_branch1_locDicaij(S_ref(:,1),S_vref(:,1),fd_array,doa_rd,f);   %%%返回瞬时定位， 和配对索引
     end
   else
      loc_shun1=-1;
   end
       
    loc_shunshi=zeros(3,1);   
    value1=size(Index1);  
    couple_array=cell(1,1);
    couple_array{1}=0;  
    if loc_shun1(1,1)~=-1
       loc_shunshi=[loc_shunshi loc_shun1];
       couple_array{1}=zeros(size(branch1_box,2),value1(2));     %%收集数据
         for i=1:value1(2)
            couple_array{1}(:,i)=branch1_box(value1(i),:)';
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
 %%%%%%%%%%%%%%%%%%%%%%分支2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 if flag_branch==2
     
     branch1_box=point_cluster{1};
     branch2_box=point_cluster{2}; 
   if size(branch1_box,2)>=6
       doa_rd=branch1_box(:,5:6);
     if flag_caiyang==1       
         tao_array=branch1_box(:,1);
         tao_array(:,2)=branch1_box(:,7);
         [loc_shun1,Index1]=Qin3_branch1_locGaocai(S_ref(:,1),tao_array,doa_rd);   %%%返回瞬时定位， 和配对索引   
     else
         fd_array=branch1_box(:,3);
         fd_array(:,2)=branch1_box(:,7);
         f=f_array;
         [loc_shun1,Index1]=Qin3_branch1_locDicaij(S_ref(:,1),S_vref(:,1),fd_array,doa_rd,f);   %%%返回瞬时定位， 和配对索引
     end
   else
      loc_shun1=-1;
   end
   
   if branch2_box(1,1) ~=-1
       if (Q-1)==1
         tao_array=branch2_box(:,1);
         tao_array(:,2)=branch2_box(:,5);
         fd_array=branch2_box(:,3);
         fd_array(:,2)=branch2_box(:,5);
       elseif (Q-1)==2
         tao_array=branch2_box(:,1);
         tao_array(:,2)=branch2_box(:,5);
         tao_array(:,3:4)=branch2_box(:,9:10);   
         
         fd_array=branch2_box(:,3);
         fd_array(:,2)=branch2_box(:,7);
         fd_array(:,3:4)=branch2_box(:,9:10);   
       elseif (Q-1)==3
         tao_array=branch2_box(:,1);
         tao_array(:,2)=branch2_box(:,5);
         tao_array(:,3)=branch2_box(:,9);  
         tao_array(:,4:6)=branch2_box(:,13:15);  
         
         fd_array=branch2_box(:,3);
         fd_array(:,2)=branch2_box(:,7);
         fd_array(:,3)=branch2_box(:,11); 
         fd_array(:,4:6)=branch2_box(:,13:15);
       end
       if flag_caiyang==1       %%%分为高采样 低采样 不同定位形式
         s_xyz=S_xyz;   %%%将所有散射体送入
         [loc_shun2,Index2]=Qin3_branch2_locGaocai(s_xyz,S_ref(:,1),tao_array,doa_ref,Q);   %%%返回瞬时定位， 和配对索引   
       else
         s_xyz=S_xyz;   %%%将所有散射体送入  除了参考散射体
         f=f_array;
         s_v=S_varray;
         [loc_shun2,Index2]=Qin3_Branch2_locDicaij(s_xyz,S_ref(:,1),s_v,S_vref(:,1),fd_array,doa_ref,f,Q);   %%%返回瞬时定位， 和配对索引
      end
   else
       loc_shun2=-1;
   end

loc_shunshi=zeros(3,1);   
value1=size(Index1); value2=size(Index2); 
couple_array=cell(1,3);
couple_array{1}=0;   couple_array{2}=0;  
if loc_shun1(1,1)~=-1
   loc_shunshi=[loc_shunshi loc_shun1];
   couple_array{1}=zeros(size(branch1_box,2),value1(2));     %%收集数据
     for i=1:value1(2)
        couple_array{1}(:,i)=branch1_box(Index1(i),:)';
     end
end
if loc_shun2(1,1)~=-1
   loc_shunshi=[loc_shunshi loc_shun2];
   if (Q-1)==2
      couple_array{2}=zeros(size(branch2_box,2),value2(2));     %%收集数据
   elseif (Q-1)==3
      couple_array{2}=zeros(size(branch2_box,2),value2(2));     %%收集数据
   end
     for i=1:value2(2)
        couple_array{2}(:,i)=branch2_box(Index2(i),:)';
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
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%分支3%%%%%%%%%%%%%%%%%%%%%
 if flag_branch==3
     
     branch1_box=point_cluster{1};
     branch2_box=point_cluster{2};
     branch3_box=point_cluster{3};
     
   if size(branch1_box,2) >=6 
       doa_rd=branch1_box(:,5:6);
     if flag_caiyang==1       %%%分为高采样 低采样 不同定位形式
         tao_array=branch1_box(:,1);
         tao_array(:,2)=branch1_box(:,7);
         [loc_shun1,Index1]=Qin3_branch1_locGaocai(S_ref(:,1),tao_array,doa_rd);   %%%返回瞬时定位， 和配对索引   
     else
         fd_array=branch1_box(:,3);
         fd_array(:,2)=branch1_box(:,7);
         f=f_array;
         [loc_shun1,Index1]=Qin3_Branch1_locDicaij(S_ref(:,1),S_vref(:,1),fd_array,doa_rd,f);   %%%返回瞬时定位， 和配对索引
     end
   else
      loc_shun1=-1;
   end
   
   if branch2_box(1,1) ~=-1
       if (Q-1)==1
         tao_array=branch2_box(:,1);
         tao_array(:,2)=branch2_box(:,5);
         
         fd_array=branch2_box(:,3);
         fd_array(:,2)=branch2_box(:,5);
       elseif (Q-1)==2
         tao_array=branch2_box(:,1);
         tao_array(:,2)=branch2_box(:,5);
         tao_array(:,3:4)=branch2_box(:,9:10);   
         
         fd_array=branch2_box(:,3);
         fd_array(:,2)=branch2_box(:,7);
         fd_array(:,3:4)=branch2_box(:,9:10);   
       elseif (Q-1)==3
         tao_array=branch2_box(:,1);
         tao_array(:,2)=branch2_box(:,5);
         tao_array(:,3)=branch2_box(:,9);  
         tao_array(:,4:6)=branch2_box(:,13:15);  
         
         fd_array=branch2_box(:,3);
         fd_array(:,2)=branch2_box(:,7);
         fd_array(:,3)=branch2_box(:,11); 
         fd_array(:,4:6)=branch2_box(:,13:15);
       end
       if flag_caiyang==1       %%%分为高采样 低采样 不同定位形式
         s_xyz=S_xyz;   %%%将所有散射体送入 除去参考散射体
         [loc_shun2,Index2]=Qin3_branch2_locGaocai(s_xyz,S_ref(:,1),tao_array,doa_ref,Q);   %%%返回瞬时定位， 和配对索引   
       else
         s_xyz=S_xyz;   %%%将所有散射体送入 除去参考散射体
         f=f_array;
         s_v=S_varray;
         [loc_shun2,Index2]=Qin3_Branch2_locDicaij(s_xyz,S_ref(:,1),s_v,S_vref(:,1),fd_array,doa_ref,f,Q);   %%%返回瞬时定位， 和配对索引
      end
   else
       loc_shun2=-1;
   end

    if branch3_box(1,1) ~=-1
        if (T-1)==1
         tao_array=branch3_box(:,1);
         tao_array(:,2)=branch3_box(:,5);
         
         fd_array=branch3_box(:,3);
         fd_array(:,2)=branch3_box(:,5);
       elseif (T-1)==2
         tao_array=branch3_box(:,1);
         tao_array(:,2)=branch3_box(:,5);
         tao_array(:,3:4)=branch3_box(:,9:10);   
         
         fd_array=branch3_box(:,3);
         fd_array(:,2)=branch3_box(:,7);
         fd_array(:,3:4)=branch3_box(:,9:10);   
       elseif (T-1)==3
         tao_array=branch3_box(:,1);
         tao_array(:,2)=branch3_box(:,5);
         tao_array(:,3)=branch3_box(:,9);  
         tao_array(:,4:6)=branch3_box(:,13:15);  
         
         fd_array=branch2_box(:,3);
         fd_array(:,2)=branch2_box(:,7);
         fd_array(:,3)=branch2_box(:,11); 
         fd_array(:,4:6)=branch2_box(:,13:15);
        end
       s_xyz=S_xyz;   %%%将所有散射体送入 除去参考散射体
       if flag_caiyang==1       %%%分为高采样 低采样 不同定位形式
         [loc_shun3,Index3]=Qin3_branch3_locGaocai(s_xyz,S_ref(:,1),tao_array,doa_ref,T);   %%%返回瞬时定位， 和配对索引   
       else
         f=f_array;
         [loc_shun3,Index3]=Qin3_Branch3_locDicaij(s_xyz,S_ref(:,1),S_vref(:,2:end),S_vref(:,1),fd_array,doa_ref,f,T);   %%%返回瞬时定位， 和配对索引
      end
   else
       loc_shun3=-1;
   end

loc_shunshi=zeros(3,1);   

couple_array=cell(1,3);
couple_array{1}=0;   couple_array{2}=0;  couple_array{3}=0;
if loc_shun1(1,1)~=-1
   loc_shunshi=[loc_shunshi loc_shun1];
   value1=size(Index1); 
   couple_array{1}=zeros(7,value1(2));     %%收集数据
     for i=1:value1(2)
        couple_array{1}(:,i)=branch1_box(Index1(i),:)';
     end
end
if loc_shun2(1,1)~=-1
   loc_shunshi=[loc_shunshi loc_shun2];
   value2=size(Index2);
   if (Q-1)==2
      couple_array{2}=zeros(size(branch2_box,2),value2(2));     %%收集数据
   else
      couple_array{2}=zeros(size(branch2_box,2),value2(2));     %%收集数据
   end
     for i=1:value2(2)
        couple_array{2}(:,i)=branch2_box(Index2(i),:)';
     end
end
if loc_shun3(1,1)~=-1
   loc_shunshi=[loc_shunshi loc_shun3];
   value3=size(Index3); 
   if (T-1)==2
      couple_array{3}=zeros(size(branch3_box,2),value3(2));     %%收集数据
   else 
      couple_array{3}=zeros(size(branch3_box,2),value3(2));     %%收集数据
   end
     for i=1:value3(2)
        couple_array{3}(:,i)=branch3_box(Index3(i),:)';
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
 

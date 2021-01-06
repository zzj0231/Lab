function [pufeng_array,Index_array]=RD_pufengsearch(RD_array,index_c,noise_power,doa_ref,si_index,flag_caiyangl,flag_noise)
%%%%情景2 经过距离谱CFAR检验后对RD谱进行谱峰搜索
%%%% RD_array 当前送入的散射体I对应的RD矩阵
%%%% index_c 若干极大值点索引号，列数据  doa_ref参考角度 ， flag_cayangl 采样标志位 1:32.5 2:16.25
%%%% noise_power 从两个CFAR中返回的噪声功率  x(k), Sn   行数据
%%%% 3:8.125  4:4.0625       0是低频： 6.199  
%%%% pufeng_array 列是数据  时差、信噪、频差、信噪、散射体索引   若等于0 证明无效删除对应散射体
%%%% Index_array  列是数据  时差序号、信噪、频差序号、信噪、参考方向(5,6) 扇区索引  主要为TDOA-DOA-FDOA数据结构体 必要时可删除

high_jingdu=[30.77e-9,61.54e-9,123.05e-9,246.15];    %%时间分辨率            已经补足4.0625
zhong=[164,83,42,22];    
di_fjindu=6.199;                                     %%低采样频率分辨率
di_tjindu=0;                                         %%低采样时间分辨率
pufeng_array=zeros(5,200);
Index_array=zeros(7,200);
p_noise=noise_power;                                 %%噪声功率   决定:采用CFAR中的噪声功率
col2=size(RD_array);

if flag_noise>=1
    tao_fenbian=high_jingdu(flag_caiyangl);
    fd_fenbian=0;                         %%需要记录高采样下高分辨率,设为0
    col=size(index_c);

  if index_c~=0
    for i=1:col(1)
        l=index_c(i);
       [C_max,k]=max(RD_array(l,:));
       %%%%计算信噪比
       if flag_caiyangl==3
           noise=p_noise(2);
       else
           noise=p_noise(1);
       end
       SNR_t=C_max/noise;
       SNR_f=0;
       %%%%对TDOA序号进行插值处理   FDOA不需要
%        if l>=2 && l<col2(1)
%            C_pre=RD_array(l-1,k);
%            C_next=RD_array(l+1,k);
%            l=l+0.5*(C_pre-C_next)/(C_pre-2*C_max+C_next);
%        end
       %%%%
       pufeng_array(1,i)=(l-zhong(flag_caiyangl))*tao_fenbian;  % *2;       % 正式运行需要将*2 剔除 因为仿真时信号产生进行了2倍频时延降为实际的一半，定位使用的是真实时延因此需要乘2
       pufeng_array(2,i)=SNR_t;             %%SNR_t已经定义
       pufeng_array(3,i)=k*fd_fenbian;   
       pufeng_array(4,i)=SNR_f;             %%SNR_f已经定义
       pufeng_array(5,i)=si_index;          %%记录散射体索引
       %%%% TDOA-DOA-FODA数据结构体所用数据  需要时再解开
       Index_array(1,i)=l;
       Index_array(2,i)=SNR_t;             
       Index_array(3,i)=k;   
       Index_array(4,i)=SNR_f;            
       Index_array(5:6,i)=doa_ref;          %%记录方向
    end
  else
      pufeng_array=-1;                    %%%若谱峰array无效 删除对应散射体
  end
else
     col=size(index_c);
   if index_c ~=0               %%%检查是否CFAR传递有效数据
     for i=1:col(1)
        k=index_c(i);
       [C_max,l]=max(RD_array(:,k));
       %%%%计算信噪比
       if flag_noise==3
           noise=p_noise(2);
       else
           noise=p_noise(1);
       end
       SNR_t=0;
       SNR_f=C_max/noise;
      %%%%对FDOA序号进行插值处理   TDOA不需要
       if k>=2 && k<col2(1)
           C_pre=RD_array(l,k-1);
           C_next=RD_array(l,k-1);
           k=k+0.5*(C_pre-C_next)/(C_pre-2*C_max+C_nest);
       end
       %%%%
       pufeng_array(1,i)=l*di_tjindu;
       pufeng_array(2,i)=SNR_t;             %%SNR_t已经定义
       pufeng_array(3,i)=k*di_fjindu;             
       pufeng_array(4,i)=SNR_f;             %%SNR_f已经定义
       pufeng_array(5,i)=si_index;          %%记录散射体索引
       %%%% 需要时再解开
       Index_array(1,i)=l;
       Index_array(2,i)=SNR_t;             %%SNR_t定义
       Index_array(3,i)=k;   
       Index_array(4,i)=SNR_f;             %%SNR_f定义
       Index_array(6:7,i)=doa_ref;          %%记录散射体索引
     end
   else
       pufeng_array=-1;
   end
   
end
if index_c(1,1)~=0
pufeng_array(:,all(pufeng_array==0,1))=[];
Index_array(:,all(Index_array==0,1))=[];
end

function point_cluster_cell=Qin3_pro_trace_couple(point_cluster,T,flag_caiyang)
 %%%应用于情景2 释放升空散射体但未检测到散射体，筛选点迹的遍历组合
 %%%point_cluster是组合定位点迹  元胞矩阵：为了应对不同矩阵数量下不同维矩阵，每个矩阵行数据组合
 %%%行是组合结果：TDOA时差、等效信噪比_tdoa、FDOA时差、等效信噪比_fdoa、方向角_ref，俯仰角_ref
 %%%T OS-CFAR门限值  S LOG-t CFAR门限值   每次杂波检验只会生成一个门限值
 %%%flag_caiyang 1: 高采样  0：低采样
 %%%point_cluster_cell 每单元矩阵 行是数据
 
cell_size=size(point_cluster{1},2);
gate_limit = 2*T-1;
if cell_size <= 7
    point_cluster_cell=cell(1,1);
    trace1=point_cluster{1};
    array_size=size(trace1);     %%%获得行数
    for i=1:array_size(1)
        if flag_caiyang==1
            SNR_T=trace1(i,2);
            if SNR_T<gate_limit 
                trace1(i,:)=zeros(1,size(trace1,2));   %%%%因为增加散射体索引因此要修改   Qin2 =7
            end                                                                      %% Qin3 至少两个时延参数才可以完成定位
        else
           SNR_f=trace1(i,4);
            if SNR_f<gate_limit
                trace1(i,:)=zeros(1,size(trace1,2));
            end  
        end
    end
    trace1(all(trace1==0,2),:)=[];
    shape=size(trace1);   row=shape(1);
    if row>0
    point_cluster_cell{1}=trace1;
    else
    point_cluster_cell{1}=-1;
    end
elseif cell_size <=10
    point_cluster_cell=cell(1,1);
    trace2=point_cluster{1};
    array2_size=size(trace2);     %%%获得行数
    %%%%处理双组组合
    for i=1:array2_size(1)
        if flag_caiyang==1
            SNR_T=[trace2(i,2),trace2(i,6)];
            SNR_ave=sum(SNR_T)/2;
            if SNR_ave<gate_limit 
                trace2(i,:)=zeros(1,size(trace2,2));
            end
        else
           SNR_f=[trace2(i,4),trace2(i,8)];
           SNR_fave=sum(SNR_f)/2;
            if SNR_fave<gate_limit
                trace2(i,:)=zeros(1,size(trace2,2));
            end  
        end
    end
    trace2(all(trace2==0,2),:)=[];
    shape=size(trace2);   row=shape(1);
    if row>0
    point_cluster_cell{1}=trace2;
    else
    point_cluster_cell{1}=-1;
    end
elseif cell_size <= 17
    point_cluster_cell=cell(1,1);
    trace3=point_cluster{1};
    array3_size=size(trace3);     %%%获得行数
    
    %%%%处理三组组合
    for i=1:array3_size(1)
        if flag_caiyang==1
            SNR_T=[trace3(i,2),trace3(i,6),trace3(i,10)];
            SNR_ave=sum(SNR_T)/3;
            if SNR_ave<gate_limit 
                trace3(i,:)=zeros(1,size(trace3,2));
            end
        else
           SNR_f=[trace3(i,4),trace3(i,8),trace3(i,12)];
           SNR_fave=sum(SNR_f)/3;
            if SNR_fave<gate_limit 
                trace3(i,:)=zeros(1,size(trace3,2));
            end  
        end
    end
    trace3(all(trace3==0,2),:)=[];
    shape=size(trace3);   row=shape(1);
    if row>0
    point_cluster_cell{1}=trace3;
    else
    point_cluster_cell{1}=-1;
    end
elseif cell_size <= 18
    point_cluster_cell=cell(1,1);
    trace4=point_cluster{1};
    array4_size=size(trace4);     %%%获得行数
    
 
    %%%%处理4站组合
    for i=1:array4_size(1)
        if flag_caiyang==1
            SNR_T=[trace4(i,2),trace4(i,6),trace4(i,10),trace4(i,14)];
            SNR_ave=sum(SNR_T)/4;
            if SNR_ave<gate_limit 
                trace4(i,:)=zeros(1,size(trace4,2));
            end
        else
           SNR_f=[trace4(i,4),trace4(i,8),trace4(i,12),trace4(i,16)];
           SNR_fave=sum(SNR_f)/4;
            if SNR_fave<gate_limit 
                trace4(i,:)=zeros(1,size(trace4,2));
            end  
        end
    end
    trace4(all(trace4==0,2),:)=[];
    shape=size(trace4);   row=shape(1);
    if row>0
    point_cluster_cell{1}=trace4;
    else
    point_cluster_cell{1}=-1;
    end
end
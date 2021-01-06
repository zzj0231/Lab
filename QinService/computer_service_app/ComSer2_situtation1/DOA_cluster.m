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
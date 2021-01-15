function location_shun=subcluster_loc(loc_array,M,ya,yb,k_up, k_down,r)
% 情景一未施放升空散射体的瞬时定位点迹计算_减法聚类晒选定为点
% loc_array组有效的站坐标,包含定位坐标xy和误差信任度函数 、俯仰角 4*M矩阵
% M为有效站位数
%location_shun 是元胞矩阵 cell：4*n xy坐标，误差信任度，俯仰角
% 注意： ya yb k_up k_down r 需要自行定义参数
% ya=5.5; yb=8.25; k_up=0.55; k_down=0.3; r=6.5;

global aalag_while_subcluster;    %确保在程序执行时该函数内的while(1)执行正常
aalag_while_subcluster=0;
xy=loc_array(1:2,:); %提取出来纯定位点坐标
shape=size(loc_array); col=shape(2);
D=zeros(1,col);             %创建各个定位点的密度指标

%%%%%%寻找第一个中心点
for i=1:col
    X=xy(:,i);
    Euler_pre=(xy-X);     
    Euler_pre2=sum(Euler_pre.^2);         %sum（矩阵）对列求和 sum（矩阵，2）对行求和
    Euler=(sqrt(Euler_pre2))./(0.5*ya)^2; %计算e的指数部分
    D(i)=sum(exp(-Euler));
end
[D_c1,index]=max(D);
sub_cen=zeros(2,100);      sub_col=1;      %创建中心记录矩阵
location_pre=cell(1,100);  loc_col=1;      %存放聚类结果
location_shun=cell(1,100);
sub_cen(:,1)=xy(:,index);   %先记录第一个中心点

    %%%%开始寻找其他聚类点
    D_ck=D_c1;
while(1)
    Euler_pre=(xy-sub_cen(:,sub_col));     
    Euler_pre2=sum(Euler_pre.^2);         %sum（矩阵）对列求和 sum（矩阵，2）对行求和
    Euler=(sqrt(Euler_pre2))./(0.5*yb)^2; %计算e的指数部分
    D=D-D_ck.*exp(-Euler);                %得到新的一轮的密度指标
    
   while(1)
       [D_ck,index]=max(D);                  %得到新的最大值和索引位
    if (D_ck<k_down*D_c1)
        flag=-1;
    elseif (D_ck>k_up*D_c1)
        sub_cen(:,sub_col+1)=xy(:,index);
        sub_col=sub_col+1;
        flag=1;
    else
        d=sub_cen-xy(:,index);
        d2=sum(d.^2);
        d_min=min(sqrt(d2));
        if (d_min/ya+D_ck/sqrt(xy(:,index)'*xy(:,index)))>=1
              sub_cen(:,sub_col+1)=xy(:,index);
              sub_col=sub_col+1;
              flag=1;
        else
            D(index)=0;
            flag=0;
            continue;
        end
    end
    if (flag==-1 || flag==1)
        break;   %跳出内循环
    end
   end
   if flag==-1
       aalag_while_subcluster=1;
       break;   %%跳出外循环
   elseif flag ==1
       flag=0;
       continue;
   end
end

for i=1:sub_col               %%%开始依据中心点构建分组
    cen=sub_cen(:,i);
    dist_pre=sum((xy-cen).^2);
    dist=sqrt(dist_pre);
    for j=1:col
        if (dist(j)<=r)       %%r需要自己定义
            location_pre{i}(:,loc_col)=loc_array(:,j);
            loc_col=loc_col+1;
        end
    end
    loc_col=1;
end
loc_col=1;                   %%location_shun继续使用

for i=1:sub_col              %%删除数据过少的矩阵
    cel=location_pre{i};
    num=size(cel);
    if num(2)>=(M*2/3)
        location_shun{loc_col}=cel;
        loc_col=loc_col+1;
    end
end
location_shun(loc_col:end)=[];     %%删除空矩阵
end
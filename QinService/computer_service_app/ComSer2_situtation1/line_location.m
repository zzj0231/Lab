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
function station_xyz=LBM_XYZ(coordinate_LBH,cen_state)
% 情景一未施放升空散射体的瞬时定位点迹计算_UCA大地坐标转全局标准站心
% coordinate_LBH各个站位经纬度坐标，数据类型：矩阵
% cen_state 参考站心 L0_B0_H0 1*3维
% 列代表各个站位经纬度,单位是弧度

shape=size(coordinate_LBH); %获取矩阵维数，shape(2)记录列数
L0=cen_state(1); B0=cen_state(2); H0=cen_state(3); 
R=6372566;  %地球半径
station_xyz=zeros(shape(1),shape(2));
for col=1:shape(2)
    LBH=coordinate_LBH(:,col);
    station_xyz(1,col)=(R+LBH(3))*cos(LBH(2))*sin(LBH(1)-L0);
    station_xyz(2,col)=(R+LBH(3))*sin(LBH(2))*cos(B0)-(R+LBH(3))*cos(LBH(2))*sin(B0)*cos(LBH(1)-L0);
    station_xyz(3,col)=(R+LBH(3))*cos(LBH(2))*cos(B0)*cos(LBH(1)-L0)+(R+LBH(3))*sin(LBH(2))*sin(B0)-R-H0;
end

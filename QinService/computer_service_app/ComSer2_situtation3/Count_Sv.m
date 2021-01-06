function S_vxyz=Count_Sv(S_xyz,vsi_array,fai_array,v_uca)
% 情景2 施放升空散射体 计算散射体瞬时速度
% S_xyz各个站位经空间坐标，数据类型：矩阵  列：数据维
% vsi_array 数据结果帧传回的各个散射体的瞬时速度
% fai_array 数据结果帧传回的各个散射体的瞬时航向角   单位是度

shape=size(S_xyz); %获取矩阵维数，shape(2)记录列数
S_vxyz=zeros(3,shape(2));
for col=1:shape(2)
    S_vxyz(1,col)=vsi_array(col)*sin(fai_array(col)*pi/180)-v_uca(1);
    S_vxyz(2,col)=vsi_array(col)*cos(fai_array(col)*pi/180)-v_uca(2);
    v2_xy=S_vxyz(1,col)^2+S_vxyz(2,col)^2;
    S_vxyz(3,col)=sqrt(vsi_array(col)^2-v2_xy);      %%需要间隔时间才可以计算出Z速度
end

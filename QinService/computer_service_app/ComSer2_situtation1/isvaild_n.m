function is=isvaild_n(station_xyz,ref_xyz)
% 情景一未施放升空散射体的瞬时定位点迹计算_第n组有效站位实时判断
% station_xyz各个站位空间坐标，数据类型：当前站位坐标
%ref_xyz 当前第n组测量下的参考站
% 有效判断原则：当前站与上一站相比若距离大于20m 则判断为有效 is=0无效 is=1 有效

    pre=station_xyz-ref_xyz;
    dis=sqrt(pre'*pre);     %
    if dis>=20
        is=1;
    else
        is=0;
    end
end
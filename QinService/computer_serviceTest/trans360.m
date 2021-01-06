function du=trans360(xyz,du)
%%%%情景一、二、三通用 用于进行角度转换
%%%%s_xyz  某个散射体的空间坐标 列数据
%%%%index_Rd 是索引号

if du<0 && xyz(1)>0
    du=du+360;
elseif du<0 && xyz(1)<0
    du=du+180;
elseif du>0 && xyz(1)<0
    du=du+180;
end
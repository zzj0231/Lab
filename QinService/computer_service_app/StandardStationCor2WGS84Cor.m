function [wgs_loc] = StandardStationCor2WGS84Cor(ssc_loc,WGS_sta)
%UNTITLED5 此处显示有关此函数的摘要
%   此处显示详细说明
% ssc_loc：定位结果，3*1
% WGS_sta：站址经纬坐标，1*3
L0 = WGS_sta(1,1);
B0 = WGS_sta(1,2);
H0 = WGS_sta(1,3);
G0 = [-sind(L0),-cosd(L0)*sind(B0),cosd(L0)*cosd(B0);
       cosd(L0),-sind(L0)*sind(B0),sind(L0)*cosd(B0);
       0,       cosd(B0),           sind(B0)];
R0 = 6372566;
X = G0*[ssc_loc(1,1);ssc_loc(2,1);ssc_loc(3,1) - H0] + [R0*cosd(B0)*cosd(L0);R0*cosd(B0)*sind(L0);R0*sind(B0)];
wgs_loc(1,1) = atand(X(2)/X(1));
wgs_loc(1,2) = atand(X(3)/sqrt(X(1)^2 + X(2)^2));
wgs_loc(1,3) = sqrt(X(1)^2 + X(2)^2)/cosd(wgs_loc(1,2)) - R0;

end


function [S_Aoa,index_Rd]=Select_Rd(s_xyz)
%%%%情景二 用于选择RD数据区域  一共12个区域
%%%%s_xyz  某个散射体的空间坐标 列数据
%%%%index_Rd 是索引号
S_Aoa=zeros(2,1);
fai_si=atan(s_xyz(2)/s_xyz(1))*180/pi;
o_si=sqrt(s_xyz'*s_xyz);
theta=asin(s_xyz(3)/o_si)*180/pi;    %%%保留fai 等于﹣的权利
if fai_si<0 && s_xyz(1)>0
    fai_si=fai_si+360;
elseif fai_si<0 && s_xyz(1)<0
    fai_si=fai_si+180;
elseif fai_si>0 && s_xyz(1)<0
    fai_si=fai_si+180;
end

S_Aoa(1,1)=fai_si;
S_Aoa(2,1)=theta;

if fai_si<180
    if fai_si<90
       if fai_si>=0 && fai_si<30
           index_Rd=1;
       elseif fai_si>=30 && fai_si<60
           index_Rd=2;
       else
           index_Rd=3;
       end
    else
       if fai_si>=90 && fai_si<120
           index_Rd=4;
       elseif fai_si>=120 && fai_si<150
           index_Rd=5;
       else
           index_Rd=6;
       end
    end
else
    if fai_si<270
        if fai_si>=180 && fai_si<210
           index_Rd=7;
       elseif fai_si>=210 && fai_si<240
           index_Rd=8;
       else
           index_Rd=9;
       end
    else
       if fai_si>=270 && fai_si<300
           index_Rd=10;
       elseif fai_si>=300 && fai_si<330
           index_Rd=11;
       else
           index_Rd=12;
       end
    end
end


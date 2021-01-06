function [fangwei_dir, yang_dir, fangwei_san, yang_san,tao_time,fd,P_dir,P_sca] = Situtation(R_xyz,T_xyz,S_xyz,T_v,f_carry)
     T_size = size(T_xyz,2);
     S_size =  size(S_xyz,2);
     fd =zeros(T_size,S_size);
     tao_time = zeros(T_size,S_size);
     radian_to_angle = 57.2958;
     P_t =20; Gt=10; sca_area=1; Gr=1;   Gdeta=10;             %%  Gt:10dB Gr:0dB gain of the tennal 
     lamda = (3e8)/(f_carry);
     P_dir = zeros(1,T_size);
     P_sca = zeros(T_size,S_size);
     for i=1:T_size              % 获得目标的直达波方向、目标直达波俯仰角
         angle_fangwei= trans360(T_xyz(:,i)-R_xyz(:,1),atan((T_xyz(2,i)-R_xyz(2,1))/(T_xyz(1,i)-R_xyz(1,1)))* radian_to_angle);
         fangwei_dir(1,i) = angle_fangwei;

         L=sqrt((T_xyz(2,i)-R_xyz(2,1)).^2+(T_xyz(1,i)-R_xyz(1,1)).^2);
         dis_T_R = sqrt(L^2+(T_xyz(3)-R_xyz(3,1))^2);
         
        P_dir(i) = (P_t*Gt*lamda^2*Gr)/(4*pi*dis_T_R)^2*1/Gdeta;         %% power of the direct wave
         
         angle_yang = atan(T_xyz(3,i)/L)*radian_to_angle;
         yang_dir(1,i) = angle_yang;
     end
     
     for i=1:S_size             % 获得散射体的方向、散射体俯仰角
         angle_fangwei= trans360(S_xyz(:,i)-R_xyz(:,1),atan((S_xyz(2,i)-R_xyz(2,1))/(S_xyz(1,i)-R_xyz(1,1)))* radian_to_angle);
         fangwei_san(1,i) = angle_fangwei;
         L=sqrt((S_xyz(2,i)-R_xyz(2,1)).^2+(S_xyz(1,i)-R_xyz(1,1)).^2);
         angle_yang = atan(S_xyz(3,i)/L)*radian_to_angle;
         yang_san(1,i) = angle_yang;
     end
     
   for i=1:T_size             %  获得多个目标直射到接收站与从散射体反射到接收站的时域；目标数量反映在时延矩阵的行上
        T= T_xyz(:,i);
       for j=1:S_size
          L=sqrt((T(2)-R_xyz(2,1)).^2+(T(1)-R_xyz(1,1)).^2);
          dis_T_R = sqrt(L^2+(T(3)-R_xyz(3,1))^2);
          
          L2=sqrt((T(2)-S_xyz(2,j)).^2+(T(1)-S_xyz(1,j)).^2);
          dis_T_S=sqrt(L2^2+(T(3)-S_xyz(3,j))^2);
          
          L3=sqrt((R_xyz(2,1)-S_xyz(2,j)).^2+(R_xyz(1,1)-S_xyz(1,j)).^2);
          dis_R_S=sqrt(L3^2+(R_xyz(3,1)-S_xyz(3,j))^2);
          
          P_sca(i,j)=(P_t*Gt)/((4*pi)*dis_T_S^2)*sca_area*1/(4*pi*(dis_R_S)^2)*lamda^2*Gr/(4*pi); 
          
          tao_time(i,j)=abs(dis_T_R-(dis_T_S+dis_R_S))/300000000;
          fd(i,j) = -1/0.2*(((T'-S_xyz(:,j)')*T_v)/dis_T_S-((T'-R_xyz')*T_v)/dis_T_R);
       end
   end
end

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
end
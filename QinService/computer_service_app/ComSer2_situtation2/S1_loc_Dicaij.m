function [Loc_xyz,Index]=S1_loc_Dicaij(S_xyz,S_varray,fd_array,doa_ref,f_array)
 %%%情景2 释放升空散射体但未检测到散射体， FDOA-DOA组合中只有一个FDOA的定位解算  目标静止
 %%%fd_array 频差矩阵 行代表组合数据，单个散射体是单列
 %%%行是组合结果：TDOA时差、等效信噪比_tdoa、FDOA时差、等效信噪比_fdoa、方向角_ref，俯仰角_ref
 %%%s_xyz： 散射体坐标: 列数据  S_varray: 列数据  f_array是各个散射体处的频率
 %%%Loc_xyz 矩阵 列是坐标数据
 
 %%%设置初始参数
 c=3e8;
 
 X1=zeros(3,1);         %%%设置初始迭代目标坐标
 dis=20;               %%%设置初始OT射线长度
 G=zeros(3,1);         %%%初始G
 fdoa_size=size(fd_array);
 Loc_xyz=zeros(3,200);
 Index=zeros(1,200);    index=0;
 
 
 if (doa_ref(1)>=0 &&doa_ref(1)<=90) || (doa_ref(1)>=270 &&doa_ref(1)<=360)
     X1(1)=dis/sqrt((1+tan(doa_ref(1)*pi/180)^2)*(1+tan(doa_ref(2)*pi/180)^2));
     X1(2)=(dis*tan(doa_ref(1)*pi/180))/sqrt((1+tan(doa_ref(1)*pi/180)^2)*(1+tan(doa_ref(2)*pi/180)^2));
 else
     X1(1)=-dis/sqrt((1+tan(doa_ref(1)*pi/180)^2)*(1+tan(doa_ref(2)*pi/180)^2));
     X1(2)=-(dis*tan(doa_ref(1)*pi/180))/sqrt((1+tan(doa_ref(1)*pi/180)^2)*(1+tan(doa_ref(2)*pi/180)^2));
 end
     X1(3)=(dis*tan(doa_ref(2)*pi/180))/sqrt((1+tan(doa_ref(2)*pi/180)^2));

     
count=0;   I=diag([1,1,1]);
for i=1:fdoa_size(1)
X=X1;  detaX=6;
s_index=fd_array(i,2);
s_xyz=S_xyz(:,s_index);
d1=sqrt(s_xyz(:,1)'*s_xyz(:,1));
f_index=fd_array(i,2);
lamda1=c/f_array(f_index);
s_varray=S_varray(:,s_index);  v=1;
  %%%%迭代位置求解%%%%
for n=1:20
    X_s1_dis=(X-s_xyz(:,1))'*(X-s_xyz(:,1));
    %X_dis=X'*X;
    X_xy=X(1:2)'*X(1:2);
    Lst=sqrt(X_s1_dis);
   % Lt=sqrt(X_dis);
    Ltxy=sqrt(X_xy);
    
    G(1)=((X-s_xyz(:,1))'*s_varray(:,1)/Lst-s_xyz(:,1)'*s_varray(:,1)/d1)/lamda1-fd_array(i,1);    %%%更新G
    G(2)=X(2)/X(1)-tan(doa_ref(1)*pi/180);
    G(3)=X(3)/sqrt(X(1:2)'*X(1:2))-tan(doa_ref(2)*pi/180);
  
    Afx1=(s_varray(1,1)*Lst^2-(X(1)-s_xyz(1,1))*(X-s_xyz(:,1))'*s_varray(:,1))/(lamda1*Lst^3);
    Afy1=(s_varray(2,1)*Lst^2-(X(2)-s_xyz(2,1))*(X-s_xyz(:,1))'*s_varray(:,1))/(lamda1*Lst^3);
    Afz1=(s_varray(3,1)*Lst^2-(X(3)-s_xyz(3,1))*(X-s_xyz(:,1))'*s_varray(:,1))/(lamda1*Lst^3);
    A=[Afx1,Afy1,Afz1;
        -X(2)/X(1)^2,1/X(1),0;
         -X(1)*X(3)/Ltxy^3,-X(2)*X(3)/Ltxy^3,1/Ltxy];                %%%更新X
    An=(A'*A);
%%%%构建阻尼系数lamda%%%
    if n<=1
         Max=max(An); ajj=Max(1);
         lamda=5*ajj;
    end
    X_pre=X; detaX_pre=detaX;                  %%记录上一个X_pre
    detaX=(An+lamda*I)\(-A'*G);
    X=X+detaX;     %%通过最小阻尼算法获得Xn+1  这里的G相当于f 即误差函数
    %%%%阻尼最小二乘的置信区域Q计算
    G_pre=G; F=0.5*(G_pre'*G_pre);
    G(1)=((X-s_xyz(:,1))'*s_varray(:,1)/Lst-s_xyz(:,1)'*s_varray(:,1)/d1)/lamda1-fd_array(i,1);    %%%更新G
    G(2)=X(2)/X(1)-tan(doa_ref(1)*pi/180);
    G(3)=X(3)/sqrt(X(1:2)'*X(1:2))-tan(doa_ref(2)*pi/180);
    F_2=0.5*(G'*G);   L_cha=0.5*detaX'*(lamda*detaX-A'*G);
    Q=(F-F_2)/L_cha;
  if(Q>0)
      lamda=max([1/3,1-(2*Q-1)^3]);
      v=2;
  else
      X=X_pre; detaX=detaX_pre;
      lamda=lamda*v;   v=2*v;
  end
    detax=sqrt(detaX'*detaX);
    if detax<=5.2
        if (doa_ref(2)>=0 && X(3)>=0) || (doa_ref(2)<0 && X(3)<0)
            count=count+1;
            Index(index+1)=i;   index=index+1;
            Loc_xyz(:,count)=X;
            break;
        end
    end
end
end
Loc_xyz(:,all(Loc_xyz==0,1))=[];     %%清除空列
if size(Loc_xyz,2)<1
    Loc_xyz = -1;
end
Index(:,all(Index==0,1))=[];     %%清除空列
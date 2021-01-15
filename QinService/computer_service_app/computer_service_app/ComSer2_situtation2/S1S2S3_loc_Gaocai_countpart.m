function [Loc_xyz,Index]=S1S2S3_loc_Gaocai_countpart(S_xyz,tao_array,doa_ref)
 %%%情景2 释放升空散射体但未检测到散射体，TDOA-DOA组合中只有3个TDOA的定位解算
 %%%tao_array 时差矩阵 行代表组合数据，n*3列 行是时差组合 
 %%%行是组合结果：TDOA时差、等效信噪比_tdoa、FDOA时差、等效信噪比_fdoa、方向角_ref，俯仰角_ref
 %%%s_xyz： 散射体坐标: 列数据
 %%%Loc_xyz 矩阵 列是坐标数据
 
 %%%设置初始参数
 c_speed=3e8;

 X1=zeros(3,1);         %%%设置初始迭代目标坐标
 dis=-10;               %%%设置初始OT射线长度
 G=zeros(3,1);         %%%初始G
 tdoa_size=size(tao_array);
 Loc_xyz=zeros(3,200);
 Index=zeros(1,200);  index=0;
 if (doa_ref(1)>=0 &&doa_ref(1)<=90) || (doa_ref(1)>=270 &&doa_ref(1)<=360)
     X1(1)=dis/sqrt((1+tan(doa_ref(1)*pi/180)^2)*(1+tan(doa_ref(2)*pi/180)^2));
     X1(2)=(dis*tan(doa_ref(1)*pi/180))/sqrt((1+tan(doa_ref(1)*pi/180)^2)*(1+tan(doa_ref(2)*pi/180)^2));
 else
     X1(1)=-dis/sqrt((1+tan(doa_ref(1)*pi/180)^2)*(1+tan(doa_ref(2)*pi/180)^2));
     X1(2)=-(dis*tan(doa_ref(1)*pi/180))/sqrt((1+tan(doa_ref(1)*pi/180)^2)*(1+tan(doa_ref(2)*pi/180)^2));
 end
     X1(3)=(dis*tan(doa_ref(2)*pi/180))/sqrt((1+tan(doa_ref(2)*pi/180)^2));
     
  %%%%%   随机初始点  %%%%% 
%    random_inital = [50*rand(),50*rand(),0*rand()];
%      N_inital=5;
%      X_inital=zeros(3,N_inital);
%      X_inital(1,:)=-500+500* rand(1,N_inital); 
%      X_inital(2,:)=-10 +500* rand(1,N_inital); 
%      X_inital(3,:)=X1(3); 
  %%%%       %%%%%
  
count=0;  I=diag([1,1,1]);
s_xyz=zeros(3,3);
clc;
for i=1:tdoa_size(1)
    s1_index=tao_array(i,4);
    s2_index=tao_array(i,5);
    s3_index=tao_array(i,6);
    s_xyz(:,1)=S_xyz(:,s1_index);s_xyz(:,2)=S_xyz(:,s2_index);s_xyz(:,3)=S_xyz(:,s3_index);
    X0 =[0,0,5]';
    d1=sqrt((X0-s_xyz(:,1))'*(X0-s_xyz(:,1)));       %% r_s1_R
    d2=sqrt((X0-s_xyz(:,2))'*(X0-s_xyz(:,2)));       %% r_s2_R
    d3=sqrt((X0-s_xyz(:,3))'*(X0-s_xyz(:,3)));  v=1; %% r_s3_R
    
      %%%% chan 算法 %%%%%%
  C=3e8;
  A=[X0-s_xyz(:,1) X0-s_xyz(:,2) X0-s_xyz(:,3)]';  A_ni = diag([1,1,1])/A;
  K1=(s_xyz(:,1)'*s_xyz(:,1));       %% r_s1_R
  K2=(s_xyz(:,2)'*s_xyz(:,2));       %% r_s2_R
  K3=(s_xyz(:,3)'*s_xyz(:,3));  v=1; %% r_s3_R
  K0=X0'*X0; 
  con1= [-d1*C*tao_array(i,1),-d2*C*tao_array(i,2),-d3*C*tao_array(i,3)];
  con2=[-d1,-d2,-d3];
  con3=[((C*tao_array(i,1))^2+K0-K1+d1^2)/2,((C*tao_array(i,2))^2+K0-K2+d2^2)/2,((C*tao_array(i,3))^2+K0-K3+d3^2)/2];
  con4=[C*tao_array(i,1),C*tao_array(i,2),C*tao_array(i,3)];
  
  p1=con3*A_ni(1,:)'; p1_=con3*A_ni(2,:)';  p1__=con3*A_ni(3,:)';
  p2=con1*A_ni(1,:)'; p2_=con1*A_ni(2,:)';  p2__=con1*A_ni(3,:)';
  q1=con4*A_ni(1,:)';  q1_=con4*A_ni(2,:)'; q1__=con4*A_ni(3,:)';
  q2=con2*A_ni(1,:)';  q2_=con2*A_ni(2,:)'; q2__=con2*A_ni(3,:)';
  
  P1=p1+p2;  P2=p1_+p2_;  P3=p1__+p2__;
  Q1=q1+q2;  Q2=q1_+q2_;  Q3=q1__+q2__;
  
  a=(Q1^2+Q2^2+Q3^2-1);  b=-2*(Q1*(X0(1)-P1)+Q2*(X0(2)-P2)+Q3*(X0(3)-P3)); c=(X0(1)-P1)^2+(X0(2)-P2)^2+(X0(3)-P3)^2;
  syms r
  root=double(solve(a*r^2+b*r+c));
  r_RT =min(real(root));
  X1(1,1)=r_RT*Q1+P1;
  X1(2,1)=r_RT*Q2+P2;
  X1(3,1)=r_RT*Q3+P3;
  %%%%%%%%%%%%%%%%%%%%%%%% chan end %%%%%%%%%%%%%%%%%%%%%
  
%     for init=1:N_inital
%         X1=X_inital(:,init);
          X=X1; detaX=6;
  %%%%迭代位置求解%%%%
    for n=1:20
        G(1)=sqrt((X-s_xyz(:,1))'*(X-s_xyz(:,1)))-sqrt(X'*X)+d1-c_speed*tao_array(i,1);    %%%更新G
        G(2)=sqrt((X-s_xyz(:,2))'*(X-s_xyz(:,2)))-sqrt(X'*X)+d2-c_speed*tao_array(i,2);
        G(3)=sqrt((X-s_xyz(:,3))'*(X-s_xyz(:,3)))-sqrt(X'*X)+d3-c_speed*tao_array(i,3);
        X_s1_dis=(X-s_xyz(:,1))'*(X-s_xyz(:,1));
        X_s2_dis=(X-s_xyz(:,2))'*(X-s_xyz(:,2));
        X_s3_dis=(X-s_xyz(:,3))'*(X-s_xyz(:,3));
        X_dis=X'*X;
        %X_xy=X(1:2)'*X(1:2);
        Lst=sqrt(X_s1_dis);
        Lst2=sqrt(X_s2_dis);
        Lst3=sqrt(X_s3_dis);
        Lt=sqrt(X_dis);           %%% 站的初始点在(0,0,0)  目标初始点禁止设置在原点
        %Ltxy=sqrt(X_xy);
        A=[(X(1)-s_xyz(1,1))/Lst-X(1)/Lt,(X(2)-s_xyz(2,1))/Lst-X(2)/Lt,(X(3)-s_xyz(3,1))/Lst-X(3)/Lt;
            (X(1)-s_xyz(1,2))/Lst2-X(1)/Lt,(X(2)-s_xyz(2,2))/Lst2-X(2)/Lt,(X(3)-s_xyz(3,2))/Lst2-X(3)/Lt;
            (X(1)-s_xyz(1,3))/Lst3-X(1)/Lt,(X(2)-s_xyz(2,3))/Lst3-X(2)/Lt,(X(3)-s_xyz(3,3))/Lst3-X(3)/Lt];                %%%更新X
        An=(A'*A);
     %%%%构建阻尼系数lamda%%%
        if n<=1
             Max=max(An); ajj=Max(1);
             lamda=5*ajj;
        end
        X_pre=X; detaX_pre=detaX;                  %%记录上一个X_pre
        detaX=(An+lamda*I)\(-A'*G);
        G_pre=G; 
      %%%% Nestrov方法 %%%%
        if n>1
            detaX= 0.2*detaX_pre+detaX;
            X_temp=X+detaX;
            G(1)=sqrt((X_temp-s_xyz(:,1))'*(X_temp-s_xyz(:,1)))-sqrt(X_temp'*X_temp)+d1-c_speed*tao_array(i,1);    %%%更新G
            G(2)=sqrt((X_temp-s_xyz(:,2))'*(X_temp-s_xyz(:,2)))-sqrt(X_temp'*X_temp)+d2-c_speed*tao_array(i,2);
            G(3)=sqrt((X_temp-s_xyz(:,3))'*(X_temp-s_xyz(:,3)))-sqrt(X_temp'*X_temp)+d3-c_speed*tao_array(i,3);
            detaX=(An+lamda*I)\(-A'*G);
            X=X+detaX;
        else
            X=X+detaX;
        end
        %%%%阻尼最小二乘的置信区域Q计算
        F=0.5*(G_pre'*G_pre);
        G(1)=sqrt((X-s_xyz(:,1))'*(X-s_xyz(:,1)))-sqrt(X'*X)+d1-c_speed*tao_array(i,1);    %%%更新G
        G(2)=sqrt((X-s_xyz(:,2))'*(X-s_xyz(:,2)))-sqrt(X'*X)+d2-c_speed*tao_array(i,2);
        G(3)=sqrt((X-s_xyz(:,3))'*(X-s_xyz(:,3)))-sqrt(X'*X)+d3-c_speed*tao_array(i,3);
        F_2=0.5*(G'*G);   L_cha=0.5*detaX'*(lamda*detaX-A'*G);
        Q=(F-F_2)/L_cha;
      if(Q>0)
          lamda=max([1/3,1-(2*Q-1)^3]);
          v=2;
      else
          X=X_pre;      detaX=detaX_pre;
          lamda=lamda*v;   v=2*v;
      end
        detax=sqrt(detaX'*detaX);
        
        %%%   画图 迭代曲线  %%%%%
            hold on,axis([-100,100,-100,100]); 
            scatter(50,50,'d','k');
            text(50+0.5,50+0.5,'T');
            scatter(X1(1),X1(2),'*','b');
            text(X1(1)-0.5,X1(2)-0.5,'Inital');
            hold on,scatter(X(1),X(2),'*'); 
            if n>1
                 plot([X_pre(1),X(1)],[X_pre(2),X(2)],'k')
                delete(text_lt);
            end
            hold on,text_lt=text([-50,X(1)-1,-50,-50],[80,X(2)-1,64,48],{['Inital: ',num2str(X1(1)),' : ',num2str(X1(2))];['',num2str(n),];['detaG: ',num2str(F_2)];['detaX: ',num2str(detax)]});
            xlabel('m / 米');
            ylabel('m / 米');
            title(' 阻尼最小二乘迭代轨迹 ')
            grid on;
            pause(0.3);
        %%%
        if F_2<=0.5
            doa_res=zeros(1,2);
            doa_res(1)=atan(X(2)/X(1))*180/pi;
            doa_res(1)=trans360(X,doa_res(1));    %%角度转换位0-360
            doa_res(2)=atan(X(3)/sqrt(X(1)^2+X(2)^2))*180/pi;
            ddoa=sqrt((doa_res(1)-doa_ref(1))*(doa_res(1)-doa_ref(1))');   %% 只凭借方向角的精度 评价数据是否正确
            if (doa_ref(2)>=0 && X(3)>=0) || (doa_ref(2)<0 && X(3)<0)
               if ddoa<=5
                count=count+1;
                Index(index+1)=i;  index=index+1;
                Loc_xyz(:,count)=X;
                break;
               end
            end
        end
    end
%     delete(text_lt);
%   end
end

Loc_xyz(:,all(Loc_xyz==0,1))=[];     %%清除空列
Index(:,all(Index==0,1))=[];     %%清除空列
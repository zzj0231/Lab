function [Loc_xyz,Index]=Qin3_branch3_locGaocai(S_xyz,S_ref,tao_array,doa_ref,Ts)
 %%%情景3 释放升空散射体检测到散射体 Q>1 T>11 分支3
 %%%tao_array 时差矩阵 行代表组合数据，n*3列 行是时差组合 
 %%%s_xyz： 散射体坐标: 列数据
 %%%Loc_xyz 矩阵 列是坐标数据
 
 %%%设置初始参数
 c=3e8;

 X1=zeros(3,1);         %%%设置初始迭代目标坐标
 dis=20;               %%%设置初始OT射线长度
 G=zeros(3,1);         %%%初始G
 tdoa_size=size(tao_array);
 Loc_xyz=zeros(3,200);
 Index=zeros(1,200);  index=0; v=1;
 
 if (doa_ref(1)>=0 &&doa_ref(1)<=90) || (doa_ref(1)>=270 &&doa_ref(1)<=360)    %如何处理初始点设置问题
     X1(1)=dis/sqrt((1+tan(doa_ref(1)*pi/180)^2)*(1+tan(doa_ref(2)*pi/180)^2));
     X1(2)=(dis*tan(doa_ref(1)*pi/180))/sqrt((1+tan(doa_ref(1)*pi/180)^2)*(1+tan(doa_ref(2)*pi/180)^2));
 else
     X1(1)=-dis/sqrt((1+tan(doa_ref(1)*pi/180)^2)*(1+tan(doa_ref(2)*pi/180)^2));
     X1(2)=-(dis*tan(doa_ref(1)*pi/180))/sqrt((1+tan(doa_ref(1)*pi/180)^2)*(1+tan(doa_ref(2)*pi/180)^2));
 end
     X1(3)=(dis*tan(doa_ref(2)*pi/180))/sqrt((1+tan(doa_ref(2)*pi/180)^2));

  if (Ts-1)>5                   % 2个散射体可测向，这种定位方式只有一个定位方程，暂时不使用
    count=0;  I=diag([1,1,1]);
    s_xyz=zeros(3,1);
    for i=1:tdoa_size(1)
    X=X1; detaX=6;
    s1_index=tao_array(i,2);    %%在情境三下冗余
    s_xyz(:,1)=S_xyz(:,s1_index);
    d1=sqrt(s_xyz(:,1)'*s_xyz(:,1));
    dref=sqrt(S_ref(:,1)'*S_ref(:,1));
     
      %%%%迭代位置求解%%%%
    for n=1:20
        G(1)=-(sqrt((X-S_ref(:,1))'*(X-S_ref(:,1)))-sqrt((X-s_xyz(:,1))'*(X-s_xyz(:,1)))+dref-d1)-c*tao_array(i,1);    %%%更新G
   
        X_s1_dis=(X-s_xyz(:,1))'*(X-s_xyz(:,1));
        X_sref_dis=(X-S_ref(:,1))'*(X-S_ref(:,1));
        Lst_ref=sqrt(X_sref_dis);
        Lst=sqrt(X_s1_dis);
   
        A=-[(X(1)-S_ref(1,1))/Lst_ref-(X(1)-s_xyz(1,1))/Lst,(X(2)-S_ref(2))/Lst_ref-(X(2)-s_xyz(2,1))/Lst,(X(3)-S_ref(3))/Lst_ref-(X(3)-s_xyz(3,1))/Lst;
            ];                %%%更新X
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
        G(1)=-(sqrt((X-S_ref(:,1))'*(X-S_ref(:,1)))-sqrt((X-s_xyz(:,1))'*(X-s_xyz(:,1)))+dref-d1)-c*tao_array(i,1);    %%%更新G
        F_2=0.5*(G'*G);   L_cha=0.5*detaX'*(lamda*detaX-A'*G);
        Q=(F-F_2)/L_cha;
      if(Q>0)
          lamda=max([1/3,1-(2*Q-1)^3]);
          v=2;
      else
          X=X_pre; detaX=detaX_pre;
          lamda=lamda*v;   v=2*v;
      end
      if F_2<=0.5
            doa_res=zeros(1,2);
            doa_res(1)=atan(X(2)/X(1))*180/pi;
            doa_res(1)=trans360(X,doa_res(1));    %%角度转换位0-360
            doa_res(2)=atan(X(3)/sqrt(X(1)^2+X(2)^2))*180/pi;
            ddoa=sqrt((doa_res-doa_ref)*(doa_res-doa_ref)');
%             if (doa_ref(2)>=0 && X(3)>=0) || (doa_ref(2)<0 && X(3)<0)
%                 if ddoa<=7
%                     count=count+1;
%                     Index(index+1)=i;  index=index+1;
%                     Loc_xyz(:,count)=X;
%                     break;
%                 end
%             end
            count=count+1;
            Index(index+1)=i;  index=index+1;
            Loc_xyz(:,count)=X;
            break;
       end
     end
    end
 end  

 if (Ts-1)==2
    count=0;  I=diag([1,1,1]);
    s_xyz=zeros(3,2);
    for i=1:tdoa_size(1)
    X=X1; detaX=6;
    s1_index=tao_array(i,3);   
    s2_index=tao_array(i,4);

    s_xyz(:,1)=S_xyz(:,s1_index);s_xyz(:,2)=S_xyz(:,s2_index);
     d1=sqrt(s_xyz(:,1)'*s_xyz(:,1));
     d2=sqrt(s_xyz(:,2)'*s_xyz(:,2));
     dref=sqrt(S_ref(:,1)'*S_ref(:,1));
     
      %%%%迭代位置求解%%%%
    for n=1:20
        G(1)=-(sqrt((X-S_ref(:,1))'*(X-S_ref(:,1)))-sqrt((X-s_xyz(:,1))'*(X-s_xyz(:,1)))+dref-d1)-c*tao_array(i,1);    %%%更新G
        G(2)=-(sqrt((X-S_ref(:,1))'*(X-S_ref(:,1)))-sqrt((X-s_xyz(:,2))'*(X-s_xyz(:,2)))+dref-d2)-c*tao_array(i,2);
    
        X_s1_dis=(X-s_xyz(:,1))'*(X-s_xyz(:,1));
        X_s2_dis=(X-s_xyz(:,2))'*(X-s_xyz(:,2));

        X_sref_dis=(X-S_ref(:,1))'*(X-S_ref(:,1));

        Lst_ref=sqrt(X_sref_dis);
        Lst=sqrt(X_s1_dis);
        Lst2=sqrt(X_s2_dis);

       
        A=-[(X(1)-S_ref(1,1))/Lst_ref-(X(1)-s_xyz(1,1))/Lst,(X(2)-S_ref(2))/Lst_ref-(X(2)-s_xyz(2,1))/Lst,(X(3)-S_ref(3))/Lst_ref-(X(3)-s_xyz(3,1))/Lst;
           (X(1)-S_ref(1,1))/Lst_ref-(X(1)-s_xyz(1,2))/Lst2,(X(2)-S_ref(2))/Lst_ref-(X(2)-s_xyz(2,2))/Lst2,(X(3)-S_ref(3))/Lst_ref-(X(3)-s_xyz(3,2))/Lst2;
            0 0 0];                %%%更新X
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
        G(1)=-(sqrt((X-S_ref(:,1))'*(X-S_ref(:,1)))-sqrt((X-s_xyz(:,1))'*(X-s_xyz(:,1)))+dref-d1)-c*tao_array(i,1);    %%%更新G
        G(2)=-(sqrt((X-S_ref(:,1))'*(X-S_ref(:,1)))-sqrt((X-s_xyz(:,2))'*(X-s_xyz(:,2)))+dref-d2)-c*tao_array(i,2);
        F_2=0.5*(G'*G);   L_cha=0.5*detaX'*(lamda*detaX-A'*G);
        Q=(F-F_2)/L_cha;
      if(Q>0)
          lamda=lamda*max([1/3,1-(2*Q-1)^3]);
          v=2;
      else
          X=X_pre; detaX=detaX_pre;
          lamda=lamda*v;   v=2*v;
      end

        %detax=sqrt(detaX'*detaX);
        if F_2<=1.2
            doa_res=zeros(1,2);
            doa_res(1)=atan(X(2)/X(1))*180/pi;
            doa_res(1)=trans360(X,doa_res(1));    %%角度转换位0-360
            doa_res(2)=atan(X(3)/sqrt(X(1)^2+X(2)^2))*180/pi;
            ddoa=sqrt((doa_res-doa_ref)*(doa_res-doa_ref)');
%             if (doa_ref(2)>=0 && X(3)>=0) || (doa_ref(2)<0 && X(3)<0)
%                if ddoa<=7
%                 count=count+1;
%                 Index(index+1)=i;  index=index+1;
%                 Loc_xyz(:,count)=X;
%                 break;
%                end
%             end
%                if ddoa<=7
%                 count=count+1;
%                 Index(index+1)=i;  index=index+1;
%                 X(3)=0;
%                 Loc_xyz(:,count)=X;
%                 break;
%                end
                count=count+1;
                Index(index+1)=i;  index=index+1;
                X(3)=0;
                Loc_xyz(:,count)=X;
                break;
        end

     end
     end
end    
     
 
if (Ts-1)==3
    count=0;  I=diag([1,1,1]);
    s_xyz=zeros(3,3);
    for i=1:tdoa_size(1)
    X=X1; detaX=6;
    s1_index=tao_array(i,4);
    s2_index=tao_array(i,5);
    s3_index=tao_array(i,6);
    s_xyz(:,1)=S_xyz(:,s1_index);s_xyz(:,2)=S_xyz(:,s2_index);s_xyz(:,3)=S_xyz(:,s3_index);
     d1=sqrt(s_xyz(:,1)'*s_xyz(:,1));
     d2=sqrt(s_xyz(:,2)'*s_xyz(:,2));
     d3=sqrt(s_xyz(:,3)'*s_xyz(:,3));
     dref=sqrt(S_ref(:,3)'*S_ref(:,3));
      %%%%迭代位置求解%%%%
    for n=1:20
        G(1)=-(sqrt((X-S_ref(:,1))'*(X-S_ref(:,1)))-sqrt((X-s_xyz(:,1))'*(X-s_xyz(:,1)))+dref-d1)-c*tao_array(i,1);    %%%更新G
        G(2)=-(sqrt((X-S_ref(:,1))'*(X-S_ref(:,1)))-sqrt((X-s_xyz(:,2))'*(X-s_xyz(:,2)))+dref-d2)-c*tao_array(i,2);
        G(3)=-(sqrt((X-S_ref(:,1))'*(X-S_ref(:,1)))-sqrt((X-s_xyz(:,3))'*(X-s_xyz(:,3)))+dref-d3)-c*tao_array(i,3);
        X_s1_dis=(X-s_xyz(:,1))'*(X-s_xyz(:,1));
        X_s2_dis=(X-s_xyz(:,2))'*(X-s_xyz(:,2));
        X_s3_dis=(X-s_xyz(:,3))'*(X-s_xyz(:,3));
        X_sref_dis=(X-S_ref(:,1))'*(X-S_ref(:,1));

        Lst_ref=sqrt(X_sref_dis);
        Lst=sqrt(X_s1_dis);
        Lst2=sqrt(X_s2_dis);
        Lst3=sqrt(X_s3_dis);
       
        A=-[(X(1)-S_ref(1,1))/Lst_ref-(X(1)-s_xyz(1,1))/Lst,(X(2)-S_ref(2))/Lst_ref-(X(2)-s_xyz(2,1))/Lst,(X(3)-S_ref(3))/Lst_ref-(X(3)-s_xyz(3,1))/Lst;
            (X(1)-S_ref(1,1))/Lst_ref-(X(1)-s_xyz(1,2))/Lst2,(X(2)-S_ref(2))/Lst_ref-(X(2)-s_xyz(2,2))/Lst2,(X(3)-S_ref(3))/Lst_ref-(X(3)-s_xyz(3,2))/Lst2;
            (X(1)-S_ref(1,1))/Lst_ref-(X(1)-s_xyz(1,3))/Lst3,(X(2)-S_ref(2))/Lst_ref-(X(2)-s_xyz(2,3))/Lst3,(X(3)-S_ref(3))/Lst_ref-(X(3)-s_xyz(3,3))/Lst3];                %%%更新X
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
        G(1)=-(sqrt((X-S_ref(:,1))'*(X-S_ref(:,1)))-sqrt((X-s_xyz(:,1))'*(X-s_xyz(:,1)))+dref-d1)-c*tao_array(i,1);    %%%更新G
        G(2)=-(sqrt((X-S_ref(:,1))'*(X-S_ref(:,1)))-sqrt((X-s_xyz(:,2))'*(X-s_xyz(:,2)))+dref-d2)-c*tao_array(i,2);
        G(3)=-(sqrt((X-S_ref(:,1))'*(X-S_ref(:,1)))-sqrt((X-s_xyz(:,3))'*(X-s_xyz(:,3)))+dref-d3)-c*tao_array(i,3);
        F_2=0.5*(G'*G);   L_cha=0.5*detaX'*(lamda*detaX-A'*G);
        Q=(F-F_2)/L_cha;
      if(Q>0)
          lamda=lamda*max([1/3,1-(2*Q-1)^3]);
          v=2;
      else
          X=X_pre; detaX=detaX_pre;
          lamda=lamda*v;   v=2*v;
      end

        %detax=sqrt(detaX'*detaX);
        if F_2<=0.5
            doa_res=zeros(1,2);
            doa_res(1)=atan(X(2)/X(1))*180/pi;
            doa_res(1)=trans360(X,doa_res(1));    %%角度转换位0-360
            doa_res(2)=atan(X(3)/sqrt(X(1)^2+X(2)^2))*180/pi;
            ddoa=sqrt((doa_res-doa_ref)*(doa_res-doa_ref)');
%             if (doa_ref(2)>=0 && X(3)>=0) || (doa_ref(2)<0 && X(3)<0)
%                 if ddoa<=7
%                     count=count+1;
%                     Index(index+1)=i;  index=index+1;
%                     Loc_xyz(:,count)=X;
%                     break;
%                 end
%             end
            count=count+1;
            Index(index+1)=i;  index=index+1;
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
function [Loc_xyz,Index_array,loc_x,loc_y]=S1_loc_Gaocai(S_xyz,tao_array,doa_ref)
 %%%�龰2 �ͷ�����ɢ���嵫δ��⵽ɢ���壬TDOA-DOA�����ֻ��һ��TDOA�Ķ�λ����
 %%%tao_array ʱ����� �д���������ݣ�����ɢ�����ǵ���
 %%%������Ͻ����TDOAʱ���Ч�����_tdoa��FDOAʱ���Ч�����_fdoa�������_ref��������_ref
 %%%S_xyz�� ɢ��������: ������   doa_ref������Ƕ�
 %%%Loc_xyz ���� ������������  Index_array��¼�������
 
 %%%���ó�ʼ����
 c=3e8;

 X1=zeros(3,1);         %%%���ó�ʼ����Ŀ������
 dis=20;               %%%���ó�ʼOT���߳���
 G=zeros(3,1);         %%%��ʼG
 tdoa_size=size(tao_array);
 Loc_xyz=zeros(3,200);
 Index_array=zeros(1,200);   index=0;
 
 if (doa_ref(1)>=0 &&doa_ref(1)<=90) || (doa_ref(1)>=270 &&doa_ref(1)<=360)
     X1(1)=dis/sqrt((1+tan(doa_ref(1)*pi/180)^2)*(1+tan(doa_ref(2)*pi/180)^2));
     X1(2)=(dis*tan(doa_ref(1)*pi/180))/sqrt((1+tan(doa_ref(1)*pi/180)^2)*(1+tan(doa_ref(2)*pi/180)^2));
 else
     X1(1)=-dis/sqrt((1+tan(doa_ref(1)*pi/180)^2)*(1+tan(doa_ref(2)*pi/180)^2));
     X1(2)=-(dis*tan(doa_ref(1)*pi/180))/sqrt((1+tan(doa_ref(1)*pi/180)^2)*(1+tan(doa_ref(2)*pi/180)^2));
 end
     X1(3)=(dis*tan(doa_ref(2)*pi/180))/sqrt((1+tan(doa_ref(2)*pi/180)^2));
%X1=[30.6,-27.45,1.456]';
     
count=0;   I=diag([1,1,1]);    %%%����������С���˵ĵ�λ����
for i=1:tdoa_size(1)
X=X1;  detaX=6;
s_index=tao_array(i,2);
s_xyz=S_xyz(:,s_index);
d0=sqrt(s_xyz(:,1)'*s_xyz(:,1));
v=1;

%%%%%%%%%%%%%%%%%%%%�������󣬻�õ㼣
loc_x=0;  x_index=0;
loc_y=0;  y_index=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%����λ�����%%%%
for n=1:20
    
    G(1)=sqrt((X-s_xyz)'*(X-s_xyz))-sqrt(X'*X)+d0-c*tao_array(i,1);    %%%����G
    G(2)=X(2)/X(1)-tan(doa_ref(1)*pi/180);
    G(3)=X(3)/sqrt(X(1:2)'*X(1:2))-tan(doa_ref(2)*pi/180);
    X_s1_dis=(X-s_xyz)'*(X-s_xyz);
    X_dis=X'*X;
    X_xy=X(1:2)'*X(1:2);
    Lst=sqrt(X_s1_dis);
    Lt=sqrt(X_dis);
    Ltxy=sqrt(X_xy);
    A=[(X(1)-s_xyz(1))/Lst-X(1)/Lt,(X(2)-s_xyz(1))/Lst-X(2)/Lt,(X(3)-s_xyz(3))/Lst-X(3)/Lt;
        -X(2)/X(1)^2,1/X(1),0;
         -X(1)*X(3)/Ltxy^3,-X(2)*X(3)/Ltxy^3,1/Ltxy];                %%%����X
    An=(A'*A);  
    %%%%��������ϵ��lamda%%%
    if n<=1
         Max=max(An); ajj=Max(1);
         lamda=5*ajj;
    end
    X_pre=X; 
    detaX_pre=detaX;                  %%��¼��һ��X_pre
    detaX=(An+lamda*I)\(-A'*G);
%     [max_step,step_index]=max(detaX(1:2));    %%%����·��  ��������·��
%     if step_index==1
%         detaX(2)=max_step/tan(doa_ref(1)*pi/180);
%     else
%         detaX(1)=max_step*tan(doa_ref(1)*pi/180);
%     end
    X=X+detaX;     %%ͨ����С�����㷨���Xn+1  �����G�൱��f ������
    %%%%������С���˵���������Q����
    G_pre=G; 
    F=0.5*(G_pre'*G_pre);
    G(1)=sqrt((X-s_xyz)'*(X-s_xyz))-sqrt(X'*X)+d0-c*tao_array(i,1);    %%%����G
    G(2)=X(2)/X(1)-tan(doa_ref(1)*pi/180);
    G(3)=X(3)/sqrt(X(1:2)'*X(1:2))-tan(doa_ref(2)*pi/180);
    F_2=0.5*(G'*G);   L_cha=0.5*detaX'*(lamda*detaX-A'*G);
    Q=(F-F_2)/L_cha;
  if(Q>0)
      lamda=lamda*max([1/3,1-(2*Q-1)^3]);
      v=2;
      %%%%%%%%%%%%%%%%%%
      loc_x(x_index+1)=X(1);  x_index=x_index+1;
      loc_y(y_index+1)=X(2);  y_index=y_index+1;
      %%%%%%%%%%%%%%%%%%%
  else
      X=X_pre; 
      detaX=detaX_pre;
      lamda=lamda*v;   v=2*v;
  end 
    %%%%%�ж��Ƿ�����
    if F_2<=1.2
        if (doa_ref(2)>=0 && X(3)>=0) || (doa_ref(2)<0 && X(3)<0)
            count=count+1;
            Index_array(index+1)=i;   index=index+1;
            Loc_xyz(:,count)=X;
            break;
        end
    end

end
end

Loc_xyz(:,all(Loc_xyz==0,1))=[];     %%�������
if size(Loc_xyz,2)<1
    Loc_xyz = -1;
end
Index_array(:,all(Index_array==0,1))=[];     %%�������
A_Rxyz = [0,0,5]';
A_Txyz = [120,80,5]';   flag_qin=3;     
if flag_qin==2
    A_Sxyz = [400,70,50;
              320,70,60]';
    A_Sv = [ 10,8,0;
            -10,8,0]';
    rmse2=zeros(size(A_Sxyz,2),100);
elseif flag_qin==3
    A_Txyz = [320,80,5]';
    A_Sxyz = [600,170,50;
              -500,170,60;
               5,-290,50]'; 
      A_Sv = [-25,8,0;
               55,18,0;
               0, 0 ,0]';
    rmse2=zeros(size(A_Sxyz,2)-1,100);
end
sim_value = zeros(2,100);
Scample_rate_4 = 4062500;         
Scample_rate_5 = 203125;          
Scample_rate_6 = 1000;
scampling=Scample_rate_5;    nf=1.38 ;    manual_snr = 1;             % 0：根据布置场景决定SNR  1：手动设置SNR
index_fdoa = zeros(2,100);
for ii=15:15
    SNR_dir=[ii]; 
    if flag_qin==2
        SNR_sca=[SNR_dir(1)-20,SNR_dir(1)-25];
    elseif flag_qin==3
        SNR_sca=[SNR_dir(1)-25,SNR_dir(1)-30,SNR_dir(1)-20];
    end
  for i=1:1
    [dir_signal_N2,Array_sig,doaStandDir,doaStandSca,SNR_dir,SNR_sca,delaypoint_sca,delaypoint_fdoa,dople_delay_sit]=produceArraySig(A_Txyz,A_Sxyz,nf,A_Sv ,scampling,SNR_dir,SNR_sca,manual_snr,flag_qin);
    if flag_qin==2
        [y_ref,z_sector_block]=producebeam_Qin2(Array_sig,doaStandDir(1,:),doaStandDir(2,:));
    elseif flag_qin==3
        [y_ref,z_sector_block]=producebeam_Qin3(Array_sig,doaStandSca(1,:),doaStandSca(2,:));
        [~, refsca_index] = max(SNR_sca);
        delaypoint_fdoa= fdoa_betsca(refsca_index,delaypoint_fdoa);
        dople_delay_sit = fdoa_betsca(refsca_index,dople_delay_sit);
    end
    [RD_cell,k,l] =produceRD(y_ref,z_sector_block,scampling);
    Order_L=cfar(RD_cell,scampling);
    orderL=Order_L-114;
    index_fdoa(:,i) = find_fdoa(-delaypoint_fdoa,orderL,114);
  end
  index_fdoa(:,all(index_fdoa==0,1))=[];
  rmse2(:,ii+1) = RMSE(dople_delay_sit,(-index_fdoa*6.2));
  plotRd(RD_cell{1,1},k,l);
end
%%

%    plotRd(RD_cell,k,l);
%    sig_length=32768;
%    plotFFT(Array_sig,scampling,sig_length)
% rmse_fdoa(:,7:31) = rmse2(:,7:31);
%  rmse_fdoaQin3(:,1:10) = rmse2(:,1:10);
  save('rmse_fdoaQin3')
% load('rmse_fdoa')
% rmse_fdoa(:,16:31)=rmse2(:,1:16);
%%
function [dir_signal_N2,Array_sig,doaStandDir,doaStandSca,SNR_dir,SNR_sca,delaypoint_sca,delaypoint_fdoa,fd]=produceArraySig(A_Txyz,A_Sxyz,nf,A_Sv,Scampling,snr_dir,snr_sca,manual_snr,flag_Qin)
    A_fc = 1500000000;
    A_Rxyz = [0,0,5]';
    [azimuth_dir_sit, pitch_dir_sit, azimuth_sca_sit, pitch_sca_sit,time_delay_sit,fd,P_dir,P_sca] = Situtation(A_Rxyz,A_Txyz,A_Sxyz,A_Sv,A_fc,flag_Qin);  %% Build situtation include:Receiver;Target;S
                                                                                                                      %%%   azimuth_dir:azimuth of direct wave at situtation
                                                                                                                      %%% pitch_dir:pitch angle of direct wave at situtation etc..
    A_dir_num = size(A_Txyz,2);                
    A_sca_num = size(A_Sxyz,2);                
    azimuth_dir = roundn(azimuth_dir_sit,-1);  % azimuth of direct wave  [0 360] included angle between the flat-xy and the height of Target
    pangle_dir = roundn(pitch_dir_sit,0);      % pitch angle of direct wave  [-90 90]
    
    doaStandDir=zeros(2,size(azimuth_dir,2));
    doaStandDir(1,:)=azimuth_dir;
    doaStandDir(2,:)=pangle_dir;
    
    azimuth_sca =roundn(azimuth_sca_sit,-1);   % azimuth of scatter wave [0 360] included angle between the flat-xy and the height of S
    pangle_sca = roundn(pitch_sca_sit,0);      % pitch angle of scatter wave [-90 90]
    doaStandSca=zeros(2,size(azimuth_sca,2));
    doaStandSca(1,:)=azimuth_sca;
    doaStandSca(2,:)=pangle_sca;
    fd_sca = fd;                               % Doppler shift
    
    Scample_rate_4 = 4062500;         % scamping rate_4: 4.0625MHZ
    Scample_rate_5 = 203125;          % scamping rate_5: 0.203125MHZ
    Scample_rate_6 = 1000;   
    if Scampling==Scample_rate_4
        dople_fd=123.98;
        Hd = filter8M_start1500k_stop1700;
    elseif Scampling==Scample_rate_5
        dople_fd=6.2; 
        Hd = filter1M_start100k_stop115;  
%         Hd = filter200_325_start100k_stop101;  
    elseif Scampling==Scample_rate_6
        dople_fd=0.03; 
        Hd = filter1000hz_start100hz_stop120hz;  
    else
        dople_fd=0;
        disp('无该采样率对应分辨率')
    end
    
    Scampling = Scample_rate_5;
    delaypoint_sca =  floor(time_delay_sit.*Scampling);   
    delaypoint_fdoa = floor(fd./dople_fd);   
    
    A_c = 300000000;                  % speed of light
    lambda = A_c/A_fc;    
    load('location.mat');             % load阵元xyz坐标值
    location = location/200*2.5;      % 阵元间距与波长比直接影响到谱峰搜索真实DOA的结果，取值不当真实DOA取值结果不如虚假峰值
    order  = 149;   
         
%   Hd = filter4M_start1500k_stop1700k;
    
    N  = 32768*3;                             % the length of generated IQ squence
    N2 = 32768;                              
    Power_noise =0.17*nf*(1e-23)*(30e6)*1*290;
    Pn = 10*log10(Power_noise)+30;
    Noise_sig = wgn(8,N2,Pn,'dBm','complex'); 
    
  if manual_snr ==0
    SNR_dir = 10*log10(P_dir/(Power_noise));
    SNR_sca = 10*log10(P_sca/(Power_noise));
  else
    SNR_dir = snr_dir;
    SNR_sca = snr_sca;
  end
    
    Array_sig = zeros(8,N2);
    sca_record_array = zeros(A_sca_num,N2); 
    decimation = 1;
    for i_d = 1:A_dir_num 
        %%% Generate Direct Wave %%%
        dir_signal_I = filter(Hd,wgn(1,N+order,SNR_dir(i_d)+Pn-3,'dBm')); 
        dir_signal_Q = filter(Hd,wgn(1,N+order,SNR_dir(i_d)+Pn-3,'dBm'));
        dir_signal_N = dir_signal_I(order+1:end)+1i*dir_signal_Q(order+1:end);   % generated IQ squence       length:32768*2 
       dir_signal_N2 = dir_signal_N(1:decimation:N2);                          % IQ squence decimated by 2  length:32768    
        
        tao_dir=zeros(1,8);                                                      % 阵元和阵心的信号间隔时延
        dir_array_sig = zeros(8,N2);           count=1;                          % direct signal wave through array
        for i = 1:2:16
            tao_dir(i) = (location(i,1)*cosd(azimuth_dir(i_d))*sind(pangle_dir(i_d))+location(i,2)*sind(azimuth_dir(i_d))*sind(pangle_dir(i_d))+location(i,3)*cosd(pangle_dir(i_d)))/A_c;
            dir_array_sig(count,:)=dir_signal_N2(1,1:N2)*exp(+1i*2*pi*A_c*tao_dir(i)/lambda);
            count = count+1;
        end  
                                %%% Generate Scatter Wave %%%   
        attenuation_sca=zeros(A_dir_num,A_sca_num);                             % 衰减因子
        sca_array_sig = zeros(8,N2);                                            % scatter signal wave through array 
        for i_s=1:A_sca_num
            sca_signal_N = zeros(1,N);                                          % current scatter wave point
            attenuation_sca(i_d,i_s)= 10^((SNR_sca(i_d,i_s)-SNR_dir(i_d))/20)*exp(-1i*2*pi*(A_fc+fd_sca(i_d,i_s))*delaypoint_sca(i_d,i_s)/Scampling);          % Computer factor: 10*log10(A^2)=SNR_sca(m)-SNR_dir(m)
            tao_sca = zeros(1,8);                                               % 散射波阵元和阵心的信号间隔时延
            sca_array_temple = zeros(8,N2); 
         
            sca_signal_N(1+delaypoint_sca(i_d,i_s):end) = dir_signal_N(1:N-delaypoint_sca(i_d,i_s));  % gengerate scatter wave through direct wave 
            sca_signal_N2=sca_signal_N(1:decimation:N2);                                     % decimate by 2 for scatter wave   两倍抽取 时延点为原来的两倍
            sca_record_array(i_s,:)=attenuation_sca(i_d,i_s)*sca_signal_N2(1,:).*exp(1i*2*pi*fd_sca(i_d,i_s)*([0:N2-1]/Scampling));      
            count = 0;
            for i = 1:2:16
               tao_sca(count+1) = (location(i,1)*cosd(azimuth_sca(i_s))*sind(pangle_sca(i_s))+location(i,2)*sind(azimuth_sca(i_s))*sind(pangle_sca(i_s))+location(i,3)*cosd(pangle_sca(i_s)))/A_c;
               sca_array_temple(count+1,:) = attenuation_sca(i_d,i_s)*sca_signal_N2(1,:).*exp(1i*2*pi*fd_sca(i_d,i_s)*([0:N2-1]/Scampling)).*exp(+1i*2*pi*A_c*tao_sca(count+1)/lambda);
               count = count+1;
            end
            sca_array_sig = sca_array_sig+sca_array_temple;
        end
        Array_sig=Array_sig+dir_array_sig+sca_array_sig;   
    end
        Array_sig = Array_sig+Noise_sig;   
end

function  [RD_cell,k,l] =produceRD(y_ref,z_sector_block,scampling)
%Generate RD spectrum
Scample_rate_1 = 32500000;        % scamping rate_1: 32.5MHZ
Scample_rate_2 = 16250000;        % scamping rate_2: 16.25MHZ
Scample_rate_3 = 8125000;         % scamping rate_3: 8.125MHZ
Scample_rate_4 = 4062500;         % scamping rate_4: 4.0625MHZ
Scample_rate_5 = 203125;          % scamping rate_5: 2.03125MHZ
if scampling ==Scample_rate_1
    l = 163; k=1;
elseif scampling ==Scample_rate_2
    l = 82; k=2;
elseif scampling ==Scample_rate_3
    l = 41; k=3;
elseif scampling ==Scample_rate_4
    l = 21; k=6;
elseif scampling ==Scample_rate_5
    l = 1;  k=113;
end
N2 = 32768;
   window_h =(0.35875*ones(1,N2)-0.48829*cos((2*pi*[1:N2])/N2)+0.1428*cos((4*pi*[1:N2])/N2)-0.01168*cos((6*pi*[1:N2])/N2));
   RD_cell=cell(1,12);
   z_sectorw=zeros(N2,2*k+1);
   y_l1 = zeros(1,N2);    y_l3 = zeros(1,N2);
   y_l1(1:N2-1)=y_ref(2:N2);
   y_l3(2:N2)=y_ref(1:N2-1);
   y=[conj(y_l1);conj(y_ref);conj(y_l3)];
   for i=1:1
    Rd_array = zeros(2*l+1,2*k+1);
      count =0;
      for j=-k:k
%           z_sectorw(:,count+1)=(z_sector_block(i,:).*window_h.*exp(1i*(2*pi*j*[1:N2])/(N2))).';
          z_sectorw(:,count+1)=z_sector_block(i,:).*window_h.*exp(1i*(2*pi*j*[1:N2])/(N2));
          cor_zy = xcorr(z_sectorw(:,count+1),y_ref);
          Rd_array(:,count+1)=cor_zy(N2-l:N2+l)';
          count =count+1;
      end
%       Rd_array =(y*z_sectorw);
      RD_cell{1,i} = Rd_array;
   end
end

function plotRd(Rd_plot,k,l)
xtrick = -l:1:l;
[Y,X] =meshgrid(-k:k,-l:l);

surf(X,Y,abs(Rd_plot),'Edgecolor','none')
xlabel('time delay / point')
ylabel('Dopple frequence / point');
zlabel('Amplitude/dB');
title('扇区1距离多普勒谱')
set(gca,'XLim',[min(xtrick),max(xtrick)])
end

function delay_betsca = fdoa_betsca(refsca_index,delaypoint_sca)
        delay_betsca=zeros(1,size(delaypoint_sca,2)-1);
        index =1;
        for i=1:size(delaypoint_sca,2)
            if i ~=refsca_index
               delay_betsca(index)=delaypoint_sca(i)-delaypoint_sca(refsca_index);
               index = index+1;
            end
        end
end
function [y_ref,z_sector_block]=producebeam_Qin2(Array_sig,azimuth_dir,pangle_dir)
   location=load('location.mat');
   location = location.location/200*2.5;
   A_dir_num=1;
   N2 = 32768; 
   A_c = 300000000;                  % speed of light
   A_fc = 1500000000;                  % carrier frequence:1.5G
   lambda = A_c/A_fc;
   
   Rxx = (1/8192)*Array_sig(:,1:8192)*Array_sig(:,1:8192)';
   [Rx_fvec,Rx_fval]=eig(Rxx);               %% Rxx:Eigenvalue Decomposition
   Rx_fval2=diag(Rx_fval)';                  %% Convert Feature value from array to vector 
   [~,I]=sort(Rx_fval2);                     %% Sort feature value and record order
   Rx_fvec2=fliplr(Rx_fvec(:,I));            %% feature vector array based the order of the value of feature value 
   Rx_Msca = Rxx+10*min(Rx_fval2)*diag(ones(1,8));  %% M矩阵用于散射波波束形成 
 
   Rx_Un = Rx_fvec2(:,A_dir_num+1:end);          %% Un信号子空间矩阵
   Rx_Block =  Rx_Un/(Rx_Un'*Rx_Un)*Rx_Un';      %% 用于扇区波束形成的波束形成
    
   vectorA_ref = zeros(8,1);     count = 0;         %% 直达波的阵列流型矢量          
   %%% Generate Reference Channel SignaL  MVDR  %%%
   for i = 1:2:16                                                         % default 1 direct wave
       tao_dir(i) = (location(i,1)*cosd(azimuth_dir(1))*sind(pangle_dir(1))+location(i,2)*sind(azimuth_dir(1))*sind(pangle_dir(1))+location(i,3)*cosd(pangle_dir(1)))/A_c;
       vectorA_ref(count+1)=exp(+1i*2*pi*A_c*tao_dir(i)/lambda);
       count = count+1;
   end
   Wref = (Rx_Msca\vectorA_ref(:,1))/(vectorA_ref(:,1)'/Rx_Msca*vectorA_ref(:,1));
   y_ref = Wref'*Array_sig;
 %%% Generate 12 Sector beam  BLOCK 
   tao_sca  = zeros(1,8);
   Sector_array=[0 30 60 90 120 150 180 210 240 270 300 330;
                10 10 10 10  10  10  10  10  10  10  10  10];
   VectorA_sector = zeros(8,12);                                    %% Array vector for every sector(扇区)
   Wsector_array_block = zeros(8,12);                               %% computer w based Block array
   z_sector_block = zeros(12,N2);
   for i_s =1:12
       count =0;
       for i=1:2:16
           tao_sca(count+1) = (location(i,1)*cosd(Sector_array(1,i_s))*sind(Sector_array(2,i_s))+location(i,2)*sind(Sector_array(1,i_s))*sind(Sector_array(2,i_s))+location(i,3)*cosd(Sector_array(2,i_s)))/A_c;
           count = count+1;
       end
       VectorA_sector(:,i_s) = exp(+1i*2*pi*A_c*tao_sca'./lambda);
       Wsector_array_block(:,i_s) = Rx_Block*VectorA_sector(:,i_s);
   end
   
   for i=1:12
       z_sector_block(i,:)=Wsector_array_block(:,i)'*Array_sig;             %% Sector beamforming based block
   end         
end

function [y_ref,z_sector_block]=producebeam_Qin3(Array_sig,azimuth_sca,pangle_sca)
   location=load('location.mat');             
   location = location.location/200*2.5;      
   A_dir_num=1;
   A_c = 300000000;                  % speed of light
   A_fc = 1500000000;                  % carrier frequence:1.5G
   lambda = A_c/A_fc;
   
   Rxx = (1/8192)*Array_sig(:,1:8192)*Array_sig(:,1:8192)';
   [Rx_fvec,Rx_fval]=eig(Rxx);               %% Rxx:Eigenvalue Decomposition
   Rx_fval2=diag(Rx_fval)';                  %% Convert Feature value from array to vector 
   [~,I]=sort(Rx_fval2);                     %% Sort feature value and record order
   Rx_fvec2= fliplr(Rx_fvec(:,I));            %% feature vector array based the order of the value of feature value 
   Rx_Msca = Rxx+10*min(Rx_fval2)*diag(ones(1,8));  %% M矩阵用于散射波波束形成 
 
   Rx_Un = Rx_fvec2(:,1+A_dir_num+1:end);    %% Un信号子空间矩阵 其余散射波+噪声
   Rx_Block =  Rx_Un/(Rx_Un'*Rx_Un)*Rx_Un';  %% 用于扇区波束形成的波束形成
   N2 = 32768;  
   vectorA_ref = zeros(8,1);             
   
   %%% Generate Reference Channel SignaL  MVDR  Computer Service situtation 3
      if A_dir_num>1
       [~,index_sca]=max(max(azimuth_sca));
      else
       [~,index_sca]=max(azimuth_sca);   
      end
       tao_sca  = zeros(1,8); count=0;
       for i=1:2:16
           tao_sca(count+1) = (location(i,1)*cosd(azimuth_sca(index_sca))*sind(pangle_sca(index_sca))+location(i,2)*sind(azimuth_sca(index_sca))*sind(pangle_sca(index_sca))+location(i,3)*cosd(pangle_sca(index_sca)))/A_c;
           count = count+1;
       end
       vectorA_ref=exp(+1i*2*pi*A_c*tao_sca'/lambda);
       Wref = (Rx_Msca\vectorA_ref(:,1))/(vectorA_ref(:,1)'/Rx_Msca*vectorA_ref(:,1));
       y_ref = Wref'*Array_sig;

 %%% Generate 12 Sector beam  BLOCK %%%
   tao_sca  = zeros(1,8);
   Sector_array=[0 30 60 90 120 150 180 210 240 270 300 330;
                10 10 10 10  10  10  10  10  10  10  10  10];
   VectorA_sector = zeros(8,12);                                    %% Array vector for every sector(扇区)
  
   Wsector_array_block = zeros(8,12);                               %% computer w based Block array
   z_sector_block = zeros(12,N2);
   for i_s =1:12
       count =0;
       for i=1:2:16
           tao_sca(count+1) = (location(i,1)*cosd(Sector_array(1,i_s))*sind(Sector_array(2,i_s))+location(i,2)*sind(Sector_array(1,i_s))*sind(Sector_array(2,i_s))+location(i,3)*cosd(Sector_array(2,i_s)))/A_c;
           count = count+1;
       end
       VectorA_sector(:,i_s) = exp(+1i*2*pi*A_c*tao_sca'./lambda);
       Wsector_array_block(:,i_s) = Rx_Block*VectorA_sector(:,i_s);
   end
   
   for i=1:12
       z_sector_block(i,:)=Wsector_array_block(:,i)'*Array_sig;             %% Sector beamforming based block
   end   
end

function [fangwei_dir, yang_dir, fangwei_san, yang_san,tao_time,fd,P_dir,P_sca] = Situtation(R_xyz,T_xyz,S_xyz,S_v,f_carry,flag_Qin)
     T_size = size(T_xyz,2);
     S_size =  size(S_xyz,2);
     fd =zeros(T_size,S_size);
     tao_time = zeros(T_size,S_size);
     radian_to_angle = 57.2958;
     P_t =20; Gt=10; sca_area=2; Gr=1;   Gdeta=10;             %%  Gt:10dB Gr:0dB gain of the tennal 
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
          fd(i,j) = -1/0.2*(((T'-S_xyz(:,j)')*S_v(:,j))/dis_T_S-((T'-R_xyz')*S_v(:,j))/dis_T_R);
       end
   end
   if flag_Qin==3
        if S_size==3
%         fd=[-215.56,206.84,60];
          fd=[-215.56,206.84,0];    
        else
              disp('Situtation函数: 情景3下，场景生成无对应散射体数量的多普了频域')
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

function rmse = RMSE(ref_value,sim_value)
   N = size(sim_value,2);
   rmse=zeros(size(ref_value,2),1);
   for i=1:size(ref_value,2)
     rmse(i) = sqrt(sum((sim_value(i,:)-ref_value(i)).^2)*1/N);
   end
end

function plotFFT(Array_sig,sig_sample,sig_length)
     %%% 绘制阵列信号 对应的功率谱图
     figure
     N2=sig_length;
     w1 = window('hamming',N2);
     fs1 = sig_sample;
     f1 = -fs1/2:fs1/N2:fs1/2;           %% fs1/N2是频谱分辨率
     S_fft=Array_sig(1,1:N2).* w1';
     plot(f1(1:N2)/1e3,20*log10(abs(fftshift(fft(S_fft)))));
     xlabel('KHZ');
     ylabel('Amplitude/dB');
     grid on
     xtrick = [- fs1/1e3/2,fs1/1e3/2];
     set(gca,'XLim',[min(xtrick),max(xtrick)],'XMinorGrid','on','XMinorTick', 'on','YMinorTick', 'on')
     set(gca,'Box','off')
     grid on
end

function Order_L=cfar(RD_cell,scampling)
    Scample_rate_1 = 32500000;        % scamping rate_1: 32.5MHZ
    Scample_rate_2 = 16250000;        % scamping rate_2: 16.25MHZ
    Scample_rate_3 = 8125000;         % scamping rate_3: 8.125MHZ
    Scample_rate_4 = 4062500;         % scamping rate_4: 4.0625MHZ
    Scample_rate_5 = 203125;          % scamping rate_5: 2.03125MHZ
    if scampling ==Scample_rate_1
        N_pf=200;  L_k = [163,1];
    elseif scampling ==Scample_rate_2
        N_pf=150;  L_k = [82,2];
    elseif scampling ==Scample_rate_3
        N_pf=70;
    elseif scampling ==Scample_rate_4
        N_pf=30;
    elseif scampling ==Scample_rate_5
        N_pf=180;  L_k = [1,113];
    end
    Rd_array=RD_cell{1,1};
    feng = zeros(1,L_k(1));
    for j=1:2*L_k(1)+1
        feng(j)=sum(abs(Rd_array(j,:)));   % 选取峰值所在列
    end
    [~,cfar_index]=max(feng);
    Rd_sector = abs(RD_cell{1,1});                                          % 取出第一扇区的数据
    flag_noise = NoisType(Rd_sector(cfar_index,:)',0.0001);     
    disp(['噪声类型：',sprintf('%d',flag_noise)])
    Rd_date = abs(Rd_array(cfar_index,:));
    Pf=1e-3; l =2*L_k(2)+1;  N=20;  k=floor(N*3/4);  Product_unit = 2;      % Rd_数据必须是行序列
    [Order_L,T_Sn,xk]=Qin_CFAR(Pf,Rd_date,l,flag_noise,N_pf,N,k,Product_unit); 
end

function fdoa_index = find_fdoa(ref_fdoa,order_L,L_K)
   fdoa_num = size(ref_fdoa,2);
   fdoa_index=zeros(fdoa_num,1);
 for ii=1:fdoa_num
   fdoa_real=ref_fdoa(ii);
  [error,index] = min(abs(fdoa_real-order_L'));
   if error<100
       fdoa_index(ii)=index-L_K;
   end
 end
end
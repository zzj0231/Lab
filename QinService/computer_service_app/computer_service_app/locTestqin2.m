A_Txyz =[320,120,5]';  
A_Sxyz = [600,65,50;
         -700,70,60]';
rmse=zeros(1,31);
sim_value=zeros(2,100);
Scample_rate_1 = 32500000;  
nf =1.5;  index=1;   manual_snr = 1;       % 0：根据布置场景决定SNR  1：手动设置SNR
for ii=40:40
 SNR_dir=[ii];
 SNR_sca=[SNR_dir(1)-20,SNR_dir(1)-25];
for i=1:1
    scampling = Scample_rate_1;
    doa_endfind2={[20.6,1],[6.2,174]};
    [Array_sig,doaStandDir,doaStandSca,SNR_dir,SNR_sca,SNR_sig,delaypoint_sca,time_delay_sit] = produceArraySig(A_Txyz,A_Sxyz,nf,SNR_dir,SNR_sca,scampling,manual_snr);
    [y_ref,z_sector_block]=producebeam_Qin2(Array_sig,doaStandDir(1,:),doaStandDir(2,:));
    [RD_cell,k,l] = produceRD(y_ref,z_sector_block,scampling);
    
    Struct_data = mixframe_qin2(doa_endfind2,A_Sxyz,RD_cell,doaStandSca(1,1:2));              % 生成混合数据帧函数

    pf=1e-3;  flag_yundong=0;
    [StructDate,UI_printInf,flag_isUsed]=Qin2_computerservice(Struct_data,pf,flag_yundong);   % 情景2定位函数
    loc_result = StructDate{11};
    if flag_isUsed>0
        loc_result = loc_result(1:2,:);
        index_loc = find_loc(A_Txyz(1:2,1),loc_result);
        if (index_loc~=0)
            sim_value(:,i)=loc_result(:,index_loc)';
        end
    end
end
sim_value(:,all(sim_value==0,1))=[];
rmse(index) = RMSE(A_Txyz(1:2,1),sim_value);
index=index+1;
end

sig_length=32678;
plotFFT(Array_sig,scampling,sig_length)

function [Array_sig,doaStandDir,doaStandSca,SNR_dir,SNR_sca,SNR_sig,delaypoint_sca,time_delay_sit] = produceArraySig(A_Txyz,A_Sxyz,nf,SNR_dir,SNR_sca,Scampling,manual_snr)
    location=load('location.mat');             
    location = location.location/200*2.5;      
    A_Rxyz = [0,0,5]';       
 
     A_Tv = [0,0,0]';
     A_fc = 1500000000;                  % carrier frequence:1.5G
    [azimuth_dir_sit, pitch_dir_sit, azimuth_sca_sit, pitch_sca_sit,time_delay_sit,fd,P_dir,P_sca] = Situtation(A_Rxyz,A_Txyz,A_Sxyz,A_Tv,A_fc);  %% Build situtation include:Receiver;Target;S
                                                                                                                      %%%   azimuth_dir:azimuth of direct wave at situtation
                                                                                                                      %%% pitch_dir:pitch angle of direct wave at situtation etc..
    A_dir_num = size(A_Txyz,2);                % the number of direct wave
    A_sca_num = size(A_Sxyz,2);                % the number of scatter wave
    azimuth_dir = roundn(azimuth_dir_sit,-1);   % azimuth of direct wave  [0 360] included angle between the flat-xy and the height of Target
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
    A_c = 300000000;                  % speed of light
    lambda = A_c/A_fc;
    
    Scample_rate_1=32500000;
    
    Scample_rate = Scampling;                                        
    delaypoint_sca =  floor(time_delay_sit.*Scample_rate);    

    order  = 149;  
    if (Scample_rate==Scample_rate_1)
%         Hd  = filter_65M_70_6_5m;
          Hd  = filter_32_5M_start15_stop15_8M;
    else
        disp('无该采样率对应滤波器')
    end
      
    N  = 32768*3;                             % the length of generated IQ squence
    N2 = 32768;                               % the length of IQ squence Decimated by 2 :2倍抽取
    
    Power_noise =nf*1.38*(1e-23)*(30e6)*1*290;
    Pn = 10*log10(Power_noise)+30;
    Noise_sig = wgn(8,N2,Pn,'dBm','complex'); % gengerate Noise
    
   if manual_snr ==0
    SNR_dir = 10*log10(P_dir/(Power_noise));
    SNR_sca = 10*log10(P_sca/(Power_noise));
  else
    P_dir = 10^(SNR_dir./10)*Power_noise;
    P_sca = 10.^(SNR_sca./10).*Power_noise;
    SNR_dir = 10*log10(P_dir/(Power_noise));
    SNR_sca = 10*log10(P_sca/(Power_noise));
    SNR_sig = 10*log10(sum(sum(P_dir)+sum(P_sca))/(Power_noise)); 
  end

    
   % Generate Array Signal: dir_signal + sca_signal + noise_sig 
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
        attenuation_sca=zeros(A_dir_num,A_sca_num);                                     % attenuation factor of scatter 衰减因子
        sca_array_sig = zeros(8,N2);                                            % scatter signal wave through array 
        for i_s=1:A_sca_num
            sca_signal_N = zeros(1,N);                                          % current scatter wave point
            attenuation_sca(i_d,i_s)= 10^((SNR_sca(i_d,i_s)-SNR_dir(i_d))/20)*exp(-1i*(A_fc+fd_sca(i_d,i_s))*delaypoint_sca(i_d,i_s)/Scample_rate);          % Computer factor: 10*log10(A^2)=SNR_sca(m)-SNR_dir(m)
            tao_sca = zeros(1,8);                                               % 散射波阵元和阵心的信号间隔时延
            sca_array_temple = zeros(8,N2); 
         
            sca_signal_N(1+delaypoint_sca(i_d,i_s):end) = dir_signal_N(1:N-delaypoint_sca(i_d,i_s));  % gengerate scatter wave through direct wave 
            sca_signal_N2=sca_signal_N(1:decimation:N2);                                     % decimate by 2 for scatter wave   两倍抽取 时延点为原来的两倍
            sca_record_array(i_s,:)=attenuation_sca(i_d,i_s)*sca_signal_N2(1,:).*exp(1i*fd_sca(i_d,i_s)*([0:N2-1]/Scample_rate));      
            count = 0;
            for i = 1:2:16
               tao_sca(count+1) = (location(i,1)*cosd(azimuth_sca(i_s))*sind(pangle_sca(i_s))+location(i,2)*sind(azimuth_sca(i_s))*sind(pangle_sca(i_s))+location(i,3)*cosd(pangle_sca(i_s)))/A_c;
               sca_array_temple(count+1,:) = attenuation_sca(i_d,i_s)*sca_signal_N2(1,:).*exp(1i*fd_sca(i_d,i_s)*([0:N2-1]/Scample_rate)).*exp(+1i*2*pi*A_c*tao_sca(count+1)/lambda);
               count = count+1;
            end
            sca_array_sig = sca_array_sig+sca_array_temple;
        end
        Array_sig=Array_sig+dir_array_sig+sca_array_sig;   
    end
        Array_sig = Array_sig+Noise_sig;  
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

function  [RD_cell,k,l] =produceRD(y_ref,z_sector_block,scampling)
%Generate RD spectrum
Scample_rate_1 = 32500000;        % scamping rate_1: 32.5MHZ
Scample_rate_2 = 16250000;        % scamping rate_2: 16.25MHZ
Scample_rate_3 = 8125000;         % scamping rate_3: 8.125MHZ
Scample_rate_4 = 4062500;         % scamping rate_4: 4.0625MHZ
Scample_rate_5 = 2.03125;          % scamping rate_5: 2.03125MHZ
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
   for i=1:12
      Rd_array = zeros(2*l+1,2*k+1);
      count =0;
      for j=-k:k
          z_sectorw=z_sector_block(i,:).*window_h.*exp(1i*(2*pi*j*[1:N2])/(N2));
          cor_zy = xcorr(z_sectorw,y_ref);
          Rd_array(:,count+1)=cor_zy(N2-l:N2+l);
          count =count+1;
      end
      RD_cell{1,i} = Rd_array;
   end
end
%  解析帧文件生成
function Struct_data = mixframe_qin2(doa_endfind2,A_Sxyz,RD_cell,azimuth_sca)
frame_decoder = load('E:\混合数据帧2解析\Mixframe2_2020_11_23_15_49_42_493_decode.mat');
frame_decoder = frame_decoder.Struct_data;
frame_decoder{20} = size(A_Sxyz,2);  % Num_Q
Num_T = 0;
frame_decoder{56} = Num_T;           % Num_T
Num_DOA =1;                                                                   % 直达波数量 + 散射波数量
frame_decoder{41} = Num_DOA;         % Num_DOA
% 2目标
% DOA = [azimuth_dir(1) pitch_dir_sit(1)                         
%        azimuth_dir(2) pitch_dir_sit(2)                                         
%        azimuth_sca(3) pitch_sca_sit(3) 
%        azimuth_sca(1) pitch_sca_sit(1) 
%        azimuth_sca(2) pitch_sca_sit(2)   
%         0     0
%         0     0 ];
% 1目标                                              % 修改1 若DOA数量发生改变
DOA = zeros(7,2);
for i=1:Num_DOA
    DOA(i,1) = doa_endfind2{i}(1,1);
    DOA(i,2) = doa_endfind2{i}(1,2);
end

T_info = [0,0                                                   
          0,0
          0,0
          0,0];
frame_decoder{1,65} = 1;   % 参考散射体的DOA序号
inital = 1;
for i=1:7
    frame_decoder{1,41+inital} =  DOA(i,1);                                  % 方向角0-360
    frame_decoder{1,41+inital+1}= DOA(i,2);                                  % 俯仰角0-90
    inital = inital + 2;
end
%3 散射体                                          % 修改2 若散射体数量发生改变
% Sxyz = [A_Sxyz(1,1),A_Sxyz(2,1),A_Sxyz(3,1),0,azimuth_sca(1);
%         A_Sxyz(1,2),A_Sxyz(2,2),A_Sxyz(3,2),0,azimuth_sca(2);
%         A_Sxyz(1,3) ,A_Sxyz(2,3),A_Sxyz(3,3),0,azimuth_sca(3);
%           0 , 0 , 0 , 0 ,0]';  
%2 散射体
Sxyz = [A_Sxyz(1,1),A_Sxyz(2,1),A_Sxyz(3,1),0,azimuth_sca(1);
       A_Sxyz(1,2),A_Sxyz(2,2),A_Sxyz(3,2),0,azimuth_sca(2);
           0 , 0 , 0, 0, 0;
          0 , 0 , 0 , 0 ,0]';  
inital =21;
for i=1:4
    frame_decoder{inital}=Sxyz(1,i);     % S_x
    frame_decoder{inital+1}=Sxyz(2,i);   % S_y
    frame_decoder{inital+2}=Sxyz(3,i);   % S_z
    frame_decoder{inital+3}=Sxyz(4,i);   % 速度
    frame_decoder{inital+4}=Sxyz(5,i);   % 方向
    inital = inital+5;
end
index=1;
for i=1:4                                  % 散射体序号、DOA序号 
    s_order=T_info(i,1);
    doa_order=T_info(i,2);
    frame_decoder{1,56+index} = s_order;                         
    frame_decoder{1,56+index+1} = doa_order;                      
    index = index+2;
end
for i=68:79
    frame_decoder{1,i} =RD_cell{1,i-67};   % 距离多普勒谱
end
Struct_data = frame_decoder;
end

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
     fangwei_san=0;
     yang_san=0;
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

function index = find_loc(ref_loc,sim_value)
   index=0;
  [error,index_order] = min(sqrt(sum((ref_loc-sim_value).^2)));
   if error<300
       index=index_order;
   end
end
function rmse = RMSE(ref_value,sim_value)
   N = size(sim_value,2);
   rmse = sqrt(sum(sum((sim_value-ref_value).^2))*1/N);
end

function plotFFT(Array_sig,sig_sample,sig_length)
     %%% 绘制阵列信号 对应的功率谱图
     N2=sig_length;
     w1 = window('hamming',N2);
     fs1 = sig_sample;
     f1 = -fs1/2:fs1/N2:fs1/2-1;           %% fs1/N2是频谱分辨率
     S_fft=Array_sig(1,1:N2).* w1';
     plot(f1/1e6,20*log10(abs(fftshift(fft(S_fft)))));
     xlabel('MHZ');
     ylabel('Amplitude/dB');
     grid on
     xtrick = - fs1/1e6/2:1:fs1/1e6/2;
     set(gca,'XLim',[min(xtrick),max(xtrick)],'XMinorGrid','on','XMinorTick', 'on','YMinorTick', 'on')
     set(gca,'Box','off')
     grid on
end
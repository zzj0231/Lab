% Gengerate IQ data and Process IQ data FOR COMPUTER SECOND SEVERCE
    A_Rxyz = [0,0,5]';         
    A_Txyz = [220,80,5]';     
%   A_Sxyz = [300,70,50;
%            -400,70,60;
%            5,-100,50]';   
    A_Sxyz = [400,70,50;    
             -520,70,60]';   
     A_Sv = [-15,18,0;
             -15,18,0]';
     A_fc = 1500000000;                  % carrier frequence:1.5G
    [azimuth_dir_sit, pitch_dir_sit, azimuth_sca_sit, pitch_sca_sit,time_delay_sit,fd,P_dir,P_sca] = Situtation(A_Rxyz,A_Txyz,A_Sxyz,A_Sv,A_fc);  %% Build situtation include:Receiver;Target;S
                                                                                                                      %%%   azimuth_dir:azimuth of direct wave at situtation
                                                                                                                      %%% pitch_dir:pitch angle of direct wave at situtation etc..
    A_dir_num = size(A_Txyz,2);                
    A_sca_num = size(A_Sxyz,2);                
    azimuth_dir = roundn(azimuth_dir_sit,-1);  % azimuth of direct wave  [0 360] included angle between the flat-xy and the height of Target
    pangle_dir = roundn(pitch_dir_sit,0);      % pitch angle of direct wave  [-90 90]
    
    azimuth_sca =roundn(azimuth_sca_sit,-1);   % azimuth of scatter wave [0 360] included angle between the flat-xy and the height of S
    pangle_sca = roundn(pitch_sca_sit,0);      % pitch angle of scatter wave [-90 90]
    fd_sca = fd;                               % Doppler shift
    
    Scample_rate_1 = 32500000;        % scamping rate_1: 32.5MHZ
    Scample_rate_2 = 16250000;        % scamping rate_2: 16.25MHZ
    Scample_rate_3 = 8125000;         % scamping rate_3: 8.125MHZ
    Scample_rate_4 = 4062500;         % scamping rate_4: 4.0625MHZ
    Scample_rate_5 = 203125;          % scamping rate_5: 0.203125MHZ
    
    Scampling = Scample_rate_5;
    delaypoint_sca =  floor(time_delay_sit.*Scampling)+1;   
    delaypoint_fdoa = floor(fd./6.2);   
    
    A_c = 300000000;                  % speed of light
    lambda = A_c/A_fc;    
    load('location.mat');             % load阵元xyz坐标值
    location = location/200*2.5;      % 阵元间距与波长比直接影响到谱峰搜索真实DOA的结果，取值不当真实DOA取值结果不如虚假峰值
    order  = 149;   
    Hd = filter1M_start100k_stop115k;           
%   Hd = filter8M_start1500k_stop1700k;
%   Hd = filter4M_start1500k_stop1700k;
    
    N  = 32768*3;                             % the length of generated IQ squence
    N2 = 32768;                              
    Power_noise =0.17*1.38*(1e-23)*(30e6)*1*290;
    Pn = 10*log10(Power_noise)+30;
    Noise_sig = wgn(8,N2,Pn,'dBm','complex'); 
    SNR_dir = 10*log10(P_dir/(Power_noise));
    SNR_sca = 10*log10(P_sca/(Power_noise)); 
     %%  Generate Array Signal: dir_signal + sca_signal + noise_sig 
    Array_sig = zeros(8,N2);
    sca_record_array = zeros(A_sca_num,N2); 
    decimation = 1;
    for i_d = 1:A_dir_num 
        %%% Generate Direct Wave %%%
        dir_signal_I = filter(Hd,wgn(1,N+order,SNR_dir(i_d)+Pn-3,'dBm')); 
        dir_signal_Q = filter(Hd,wgn(1,N+order,SNR_dir(i_d)+Pn-3,'dBm'));
        Noise_sig_filter =filter(Hd,Noise_sig); 
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
            attenuation_sca(i_d,i_s)= 10^((SNR_sca(i_d,i_s)-SNR_dir(i_d))/20)*exp(-1i*(A_fc+fd_sca(i_d,i_s))*delaypoint_sca(i_d,i_s)/Scampling);          % Computer factor: 10*log10(A^2)=SNR_sca(m)-SNR_dir(m)
            tao_sca = zeros(1,8);                                               % 散射波阵元和阵心的信号间隔时延
            sca_array_temple = zeros(8,N2); 
         
            sca_signal_N(1+delaypoint_sca(i_d,i_s):end) = dir_signal_N(1:N-delaypoint_sca(i_d,i_s));  % gengerate scatter wave through direct wave 
            sca_signal_N2=sca_signal_N(1:decimation:N2);                                     % decimate by 2 for scatter wave   两倍抽取 时延点为原来的两倍
            sca_record_array(i_s,:)=attenuation_sca(i_d,i_s)*sca_signal_N2(1,:).*exp(1i*fd_sca(i_d,i_s)*([0:N2-1]/Scampling));      
            count = 0;
            for i = 1:2:16
               tao_sca(count+1) = (location(i,1)*cosd(azimuth_sca(i_s))*sind(pangle_sca(i_s))+location(i,2)*sind(azimuth_sca(i_s))*sind(pangle_sca(i_s))+location(i,3)*cosd(pangle_sca(i_s)))/A_c;
               sca_array_temple(count+1,:) = attenuation_sca(i_d,i_s)*sca_signal_N2(1,:).*exp(1i*fd_sca(i_d,i_s)*([0:N2-1]/Scampling)).*exp(+1i*2*pi*A_c*tao_sca(count+1)/lambda);
               count = count+1;
            end
            sca_array_sig = sca_array_sig+sca_array_temple;
        end
        Array_sig=Array_sig+dir_array_sig+sca_array_sig;   
    end
        Array_sig = Array_sig+Noise_sig;    

%%  仿真图
     xaxi = linspace(-N2+1,N2-1,2*N2-1);   
     %%% 绘制阵列信号 对应的功率谱图
     figure
     w1 = window('hamming',N2);
     fs1 = 1000000;
     f1 = -fs1/2:fs1/N2:fs1/2-1;           %% fs1/N2是频谱分辨率
     S_fft=Array_sig(1,1:N2).* w1';
     plot(f1/1e3,20*log10(abs(fftshift(fft(S_fft)))));
     xlabel('MHZ');
     ylabel('Amplitude/dB');
     title ('8通道阵列接收设备 载波频率:1.5GHZ 带宽:3M ')
     grid on
     xtrick = - 500:1:500;
     set(gca,'XLim',[min(xtrick),max(xtrick)],'XMinorGrid','on','XMinorTick', 'on','YMinorTick', 'on')
     set(gca,'Box','off')
     grid on
 %% 参考通道波束形成 与 扇区波束形成 
   Rxx = (1/8192)*Array_sig(:,1:8192)*Array_sig(:,1:8192)';
   [Rx_fvec,Rx_fval]=eig(Rxx);               %% Rxx:Eigenvalue Decomposition
   Rx_fval2=diag(Rx_fval)';                  %% Convert Feature value from array to vector 
   [~,I]=sort(Rx_fval2);                     %% Sort feature value and record order
   Rx_fvec2=fliplr(Rx_fvec(:,I));            %% feature vector array based the order of the value of feature value 
   Rx_Msca = Rxx+10*min(Rx_fval2)*diag(ones(1,8));  %% M矩阵用于散射波波束形成 
 
   Rx_Us = Rx_fvec2(:,1:A_dir_num);              %% Us信号子空间矩阵
   Rx_Un = Rx_fvec2(:,A_dir_num+1:end);          %% Un信号子空间矩阵
   Rx_Block =  Rx_Un/(Rx_Un'*Rx_Un)*Rx_Un';      %% 用于扇区波束形成的波束形成
    
   vectorA_ref = zeros(8,1);     count = 0;         %% 直达波的阵列流型矢量          
   
   %%% Generate Reference Channel SignaL  MVDR  %%%
   % Computer Service situtation 2
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
   
   for i=1:12             %% Sector beamforming based MVDR
       z_sector_block(i,:)=Wsector_array_block(:,i)'*Array_sig;             %% Sector beamforming based block
   end         

%% Generate RD spectrum
   if Scampling ==Scample_rate_1
       l = 163; k=1;
   elseif Scampling ==Scample_rate_2
       l = 82; k=2;
   elseif Scampling ==Scample_rate_3
       l = 41; k=3;
   elseif Scampling ==Scample_rate_4
       l = 21; k=6;
   elseif Scampling ==Scample_rate_5
       l = 1;  k=113;
   end
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
%%
figure
xtrick = -l:1:l;
[Y,X] =meshgrid(-k:k,-l:l);
Rd_plot = RD_cell{1,10};
surf(X,Y,abs(Rd_plot),'Edgecolor','none')
xlabel('time delay / point')
ylabel('Dopple frequence / point');
zlabel('Amplitude/dB');
title('扇区1距离多普勒谱')
set(gca,'XLim',[min(xtrick),max(xtrick)])
%% CFAR检验
if Scampling ==Scample_rate_1
    N_pf=200;  L_k = [163,1];
elseif Scampling ==Scample_rate_2
    N_pf=150;  L_k = [82,2];
elseif Scampling ==Scample_rate_3
    N_pf=70;
elseif Scampling ==Scample_rate_4
    N_pf=30;
elseif Scampling ==Scample_rate_5
    N_pf=160;  L_k = [1,113];
end
Rd_array=RD_cell{1,7};
feng = zeros(1,L_k(2));
for j=1:2*L_k(1)+1
    feng(j)=sum(abs(Rd_array(j,:)));   % 选取峰值所在列
end
[~,cfar_index]=max(feng);
Rd_sector = abs(RD_cell{1,1});                                          % 取出第一扇区的数据
flag_noise = NoisType(Rd_sector(cfar_index,:),0.0001);     
disp(['噪声类型：',sprintf('%d',flag_noise)])
Rd_date = abs(Rd_array(cfar_index,:));
Pf=10e-4; l =2*L_k(2)+1;  N=20;  k=floor(N*3/4);  Product_unit = 2;        % Rd_数据必须是行序列
[Order_L,T_Sn,xk]=Qin_CFAR(Pf,Rd_date,l,flag_noise,N_pf,N,k,Product_unit); % l+delay_point 应该是CFAR输出的对应delya_point点
disp(['CFAR门限值：',sprintf('T:%0.2f, Sn:%0.2f',T_Sn(1),T_Sn(2))])

xtrick = -L_k(2):1:L_k(2);

value =zeros(1,size(Order_L,1));
for i=1:size(Order_L,1)
    if Order_L(i)>0
        value(i)=1;
    end
end
plot(xtrick,value);
 xlabel('time delay/ point');
 ylabel('检测结果');
title('扇区7距离谱数据恒虚警检测 杂波类型：瑞利分布')
set(gca,'XLim',[min(xtrick),max(xtrick)],'YLim',[0,1.5],'XMinorGrid','on','XMinorTick','on');
grid on
%%  解析帧文件生成
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
frame_order = 7;                                    % 修改3 若保存序号发生改变        
str = sprintf('E:\\混合数据帧2解析\\验证报告生成帧\\Mixframe2_FormQin2-S%d-T%d-序号%d_decode',A_sca_num,A_dir_num,frame_order);
save(str,'Struct_data')
%%
   corrcoef(dir_signal_N2(1,:),    y_ref)               %%计算协方差系数
   corrcoef(sca_record_array(1,:), y_ref)               %%计算协方差系数
%  corrcoef(sca_record_array(2,:), z1_2)                %%计算协方差系数
   figure
   xtrick = -35:5:35;
   subplot(2,1,1)
     w1 = window('hamming',N2);
     fs1 = 65000000;
     f1 = -fs1/2:fs1/N2:fs1/2-1;           %%频谱分辨率
     S_fft=dir_signal_N2(1,1:N2).* w1';
     plot(f1/1e6,20*log10(abs(fftshift(fft(S_fft)))));
     xlabel('MHZ');
     ylabel('Amplitude/dB');
     title ('8通道阵列接收设备 直达波信号频谱图 带宽:30MHZ ')
     set(gca,'XLim',[min(xtrick),max(xtrick)],'XMinorGrid','on','XMinorTick', 'on','YMinorTick', 'on')
     set(gca,'Box','off')
     grid on
   subplot(2,1,2)
       w1 = window('hamming',N2);
       fs1 = 65000000;
       f1 = -fs1/2:fs1/N2:fs1/2-1;           %%频谱分辨率
       S_fft=y_ref(1,1:N2).* w1';
       plot(f1/1e6,20*log10(abs(fftshift(fft(S_fft)))));
     xlabel('MHZ');
     ylabel('Amplitude/dB');
     title ('8通道阵列接收设备 直达波信号频谱图 带宽:30MHZ ')
     set(gca,'XLim',[min(xtrick),max(xtrick)],'XMinorGrid','on','XMinorTick', 'on','YMinorTick', 'on')
     set(gca,'Box','off')
     grid on
   figure
    subplot(3,1,1)
     cor_S_1 =20*log10(abs(xcorr(Array_sig(1,1:N2),y_ref)));
     plot(xaxi,cor_S_1)
     title('测试： 阵列信号与参考通道信号互相关')
%      set(gca, 'GridAlpha', 1,'GridColor',[0 0 0],'MinorGridColor',[0 0 0],'XMinorTick','on','YMinorTick','on');  % 设置透明度
     xlabel('time delay / point')
     ylabel('Amplitude/dB');
     grid on
     
    subplot(2,1,1)
     cor_S_2 =20*log10(abs(xcorr(y_ref,z_sector_block(1,:))));
     plot(xaxi,cor_S_2)
     title('参考通道信号与扇区1波束形成信号互相关')
%      set(gca, 'GridAlpha', 1,'GridColor',[0 0 0],'MinorGridColor',[0 0 0],'XMinorTick','on','YMinorTick','on');  % 设置透明度
     xlabel('time delay / point')
     ylabel('Amplitude/dB');
     set(gca,'XMinorGrid','on','XMinorTick', 'on','YMinorTick', 'on')
     set(gca,'Box','off')
     grid on
    subplot(2,1,2)
     cor_S_3 =20*log10(abs(xcorr(y_ref,z_sector_block(7,:))));
     plot(xaxi,cor_S_3)
     title('参考通道信号与扇区7波束形成信号互相关')
%      set(gca, 'GridAlpha', 1,'GridColor',[0 0 0],'MinorGridColor',[0 0 0],'XMinorTick','on','YMinorTick','on');  % 设置透明度
     xlabel('time delay / point')
     ylabel('Amplitude/dB');
     set(gca,'XMinorGrid','on','XMinorTick', 'on','YMinorTick', 'on')
     set(gca,'Box','off')
     grid on
     
  figure
     cor_S_3 =20*log10(abs(xcorr(y_ref,z_sector_block(10,:))));
     plot(xaxi,cor_S_3)
     title('测试： 参考通道信号与波束形成：扇区信号10互相关')
%      set(gca, 'GridAlpha', 1,'GridColor',[0 0 0],'MinorGridColor',[0 0 0],'XMinorTick','on','YMinorTick','on');  % 设置透明度
     xlabel('time delay / point')
     ylabel('Amplitude/dB');
     grid on
%%
   similar_sectorwave = zeros(1,12);
   for i=1:A_sca_num
       for j=1:12
           cor = corrcoef(sca_record_array(i,:), z_sector_block(j,:));                    %%计算协方差系数
           similar_sectorwave(i,j)=abs(cor(1,2));
       end
   end
   figure
   subplot(2,1,1)
     w1 = window('hamming',N2);
     fs1 = 16250000;
     f1 = -fs1/2:fs1/N2:fs1/2-1;           %%频谱分辨率
     S_fft=sca_record_array(1,1:N2).* w1';
     plot(f1,20*log10(abs(fftshift(fft(S_fft)))));
     xlabel('HZ');
     ylabel('Amplitude/dB');
     title ('8通道阵列接收设备 方位角13°散射波信号频谱图 带宽:10MHZ ')
     grid on
   subplot(2,1,2)
       w1 = window('hamming',N2);
       fs1 = 16250000;
       f1 = -fs1/2:fs1/N2:fs1/2-1;         %%频谱分辨率
       S_fft=z_sector_block(5,1:N2).* w1';
       plot(f1,20*log10(abs(fftshift(fft(S_fft)))));
       xlabel('HZ');
       ylabel('Amplitude/dB');
       title ('8通道阵列接收设备 0°扇区1波束形成信号频谱图 带宽:10MHZ ')
       grid on  
   figure
   subplot(2,1,1)
     w1 = window('hamming',N2);
     fs1 = 16250000;
     f1 = -fs1/2:fs1/N2:fs1/2-1;           %%频谱分辨率
     S_fft=sca_record_array(2,1:N2).* w1';
     plot(f1,20*log10(abs(fftshift(fft(S_fft)))));
     xlabel('HZ');
     ylabel('Amplitude/dB');
     title ('8通道阵列接收设备  方位角170°散射波信号频谱图 带宽:10MHZ ')
     grid on
   subplot(2,1,2)
       w1 = window('hamming',N2);
       fs1 = 16250000;
       f1 = -fs1/2:fs1/N2:fs1/2-1;           %%频谱分辨率
       S_fft=z_sector_block(7,1:N2).* w1';
       plot(f1,20*log10(abs(fftshift(fft(S_fft)))));
       xlabel('HZ');
       ylabel('Amplitude/dB');
       title ('8通道阵列接收设备 扇区7波束形成信号频谱图 带宽:10MHZ ')
       grid on  
   figure
   subplot(2,1,1)
     w1 = window('hamming',N2);
     fs1 = 16250000;
     f1 = -fs1/2:fs1/N2:fs1/2-1;           %%频谱分辨率
     S_fft=sca_record_array(3,1:N2).* w1';
     plot(f1,20*log10(abs(fftshift(fft(S_fft)))));
     xlabel('HZ');
     ylabel('Amplitude/dB');
     title ('8通道阵列接收设备 方位角270°直达波信号频谱图 带宽:10MHZ ')
     grid on
   subplot(2,1,2)
       w1 = window('hamming',N2);
       fs1 = 16250000;
       f1 = -fs1/2:fs1/N2:fs1/2-1;           %%频谱分辨率
       S_fft=z_sector_block(10,1:N2).* w1';
       plot(f1,20*log10(abs(fftshift(fft(S_fft)))));
       xlabel('HZ');
       ylabel('Amplitude/dB');
       title ('8通道阵列接收设备 扇区10参考通道信号频谱图 带宽:10MHZ ')
       grid on  
%%
function [fangwei_dir, yang_dir, fangwei_san, yang_san,tao_time,fd,P_dir,P_sca] = Situtation(R_xyz,T_xyz,S_xyz,S_v,f_carry)
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
          fd(i,j) = -1/0.2*(((T'-S_xyz(:,j)')*S_v(:,i))/dis_T_S-((T'-R_xyz')*S_v(:,i))/dis_T_R);
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

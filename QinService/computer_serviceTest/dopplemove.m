A_fc = 1500000000;
A_Txyz =[220,120,5]';  
A_Rxyz = [0,0,5]';
A_Sxyz = [600,65,50
         -520,70,60]';
  A_Sv = [10,8,0;
          0,0,0;]';
   nf = 1.38;
[azimuth_dir_sit, pitch_dir_sit, azimuth_sca_sit, pitch_sca_sit,time_delay_sit,fd,P_dir,P_sca] = Situtation(A_Rxyz,A_Txyz,A_Sxyz,A_Sv,A_fc);  
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
fd_sca = fd;  %

Scample_rate_6 = 1000;
Scampling = Scample_rate_6;
fd_unit = Scampling/5000;
delayfre_point=floor(fd./fd_unit);

Hd = filter1000hz_start100hz_stop120hz; 
delaypoint_sca =  floor(time_delay_sit.*Scampling);  

A_c = 300000000;                  % speed of light
lambda = A_c/A_fc;
load('location.mat');             % load阵元xyz坐标值
location = location/200*2.5;      % 阵元间距与波长比直接影响到谱峰搜索真实DOA的结果，取值不当真实DOA取值结果不如虚假峰值
order  = 149;

N  = 5000*3;                             % the length of generated IQ squence
N2 = 5000;
Power_noise =0.17*nf*(1e-23)*(30e6)*1*290;
Pn = 10*log10(Power_noise)+30;
Noise_sig = wgn(8,N2,Pn,'dBm','complex');

SNR_dir = 10*log10(P_dir/(Power_noise));
SNR_sca = 10*log10(P_sca/(Power_noise));
Array_sig = zeros(8,N2);
sca_record_array = zeros(A_sca_num,N2);
decimation = 1;
for i_d = 1:A_dir_num
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
    sca_sig=zeros(A_sca_num,N2);
    sca_array_sig = zeros(8,N2);                                           
    
    for i_s=1:A_sca_num
        sca_signal_N = zeros(1,N);                                          % current scatter wave point
        attenuation_sca(i_d,i_s)= 10^((SNR_sca(i_d,i_s)-SNR_dir(i_d))/20)*exp(-1i*2*pi*(A_fc+delayfre_point(i_d,i_s)*fd_unit)*delaypoint_sca(i_d,i_s)/Scampling);          % Computer factor: 10*log10(A^2)=SNR_sca(m)-SNR_dir(m)
        tao_sca = zeros(1,8);                                               % 散射波阵元和阵心的信号间隔时延
        sca_array_temple = zeros(8,N2);
        
        sca_signal_N(1+delaypoint_sca(i_d,i_s):end) = dir_signal_N(1:N-delaypoint_sca(i_d,i_s));  % gengerate scatter wave through direct wave
        sca_signal_N2=sca_signal_N(1:decimation:N2);                                     % decimate by 2 for scatter wave   两倍抽取 时延点为原来的两倍
        sca_sig(i_s,:) =180*attenuation_sca(i_d,i_s)*dir_signal_N2.*exp(1i*2*pi*delayfre_point(i_d,i_s)*fd_unit*([0:N2-1]/Scampling));
        sca_record_array(i_s,:)=attenuation_sca(i_d,i_s)*sca_signal_N2(1,:).*exp(1i*fd_sca(i_d,i_s)*([0:N2-1]/Scampling));
        count = 0;
        for i = 1:2:16
            tao_sca(count+1) = (location(i,1)*cosd(azimuth_sca(i_s))*sind(pangle_sca(i_s))+location(i,2)*sind(azimuth_sca(i_s))*sind(pangle_sca(i_s))+location(i,3)*cosd(pangle_sca(i_s)))/A_c;
            sca_array_temple(count+1,:) = sca_sig(i_s,:).*exp(+1i*2*pi*A_c*tao_sca(count+1)/lambda);
            count = count+1;
        end
        sca_array_sig = sca_array_sig+sca_array_temple;
    end
    Array_sig=Array_sig+dir_array_sig+sca_array_sig;
end
Array_sig = Array_sig+Noise_sig;

sig_length=5000;
plotFFT(dir_array_sig,sca_array_sig,Scampling,sig_length);
plotFFT(dir_array_sig,Array_sig,Scampling,sig_length);
sig_length2=5000;
plotDople(dir_signal_N2,sca_sig,sig_length2);


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
          fd(i,j) = -1/0.2*(((T'-S_xyz(:,j)')*S_v(:,j))/dis_T_S-((T'-R_xyz')*S_v(:,j))/dis_T_R);
       end
   end
end

function plotFFT(Array_sig1,Array_sig2,sig_sample,sig_length)
     %%% 绘制阵列信号 对应的功率谱图
     figure
     N2=sig_length;
     w1 = window('hamming',N2);
     fs1 = sig_sample;
     f1 = -fs1/2:fs1/N2:fs1/2;           %% fs1/N2是频谱分辨率
     S_fft=Array_sig1(1,1:N2).*  w1';
     S_fft2=Array_sig2(1,1:N2).* w1';
     plot(f1(1:N2)/1e3,20*log10(abs(fftshift(fft(S_fft))))+10);hold on
     plot(f1(1:N2)/1e3,20*log10(abs(fftshift(fft(S_fft2)))));hold on
     xlabel('KHZ');
     ylabel('Amplitude/dB');
     grid on
     xtrick = [- fs1/1e3/2,fs1/1e3/2];
     set(gca,'XLim',[min(xtrick),max(xtrick)],'XMinorGrid','on','XMinorTick', 'on','YMinorTick', 'on')
     set(gca,'Box','off')
     legend('原始信号','多普勒频移信号');
     grid on
end

function plotDople(sig_dir,sig_sca,sig_length)
     N2=sig_length;
     k=800;
     window_h =(0.35875*ones(1,N2)-0.48829*cos((2*pi*[1:N2])/N2)+0.1428*cos((4*pi*[1:N2])/N2)-0.01168*cos((6*pi*[1:N2])/N2));
     dopler=zeros(1,2*k+1);
     count=1;
     sig_sca=sum(sig_sca);
     for j=-k:k
         dopler(count)=sum(conj(sig_dir(1,1:N2)).*sig_sca(1,1:N2).*window_h.*exp(1i*(2*pi*j*[1:N2])/(N2)));
         count =count +1;
     end
     figure;
     plot(-k:k,20*log10(abs(dopler)));
     title('多普勒谱');
     grid on;
     figure
     xaxis =linspace(-N2,N2,sig_length*2-1);
     sig=xcorr(sig_dir(1,1:N2),sig_sca(1,1:N2));
     plot(xaxis,20*log10(abs(sig)));
     title('距离谱');
end
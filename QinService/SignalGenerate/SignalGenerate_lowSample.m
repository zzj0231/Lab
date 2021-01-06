% Gengerate IQ data and Process IQ data FOR COMPUTER SECOND SEVERCE
    A_Rxyz = [0,0,5]';       % location of Receiver column 
    A_Txyz = [250,50,5]';     % location of Target   column
    A_Sxyz = [500,70,60;
             -600,70,50;
              5,-200,50]';   % location of S_array  column
     A_Tv = [0,0,0]';
     A_fc = 1500000000;                  % carrier frequence:1.5G
    [azimuth_dir_sit, pitch_dir_sit, azimuth_sca_sit, pitch_sca_sit,time_delay_sit,fd,P_dir,P_sca] = Situtation(A_Rxyz,A_Txyz,A_Sxyz,A_Tv,A_fc);  %% Build situtation include:Receiver;Target;S
                                                                                                                      %%%   azimuth_dir:azimuth of direct wave at situtation
                                                                                                                      %%% pitch_dir:pitch angle of direct wave at situtation etc..
    A_dir_num = size(A_Txyz,2);                % the number of direct wave
    A_sca_num = size(A_Sxyz,2);                % the number of scatter wave
   azimuth_dir = roundn(azimuth_dir_sit,-1);   % azimuth of direct wave  [0 360] included angle between the flat-xy and the height of Target
    pangle_dir = roundn(pitch_dir_sit,0);      % pitch angle of direct wave  [-90 90]
    
    azimuth_sca =roundn(azimuth_sca_sit,-1);   % azimuth of scatter wave [0 360] included angle between the flat-xy and the height of S
    pangle_sca = roundn(pitch_sca_sit,0);      % pitch angle of scatter wave [-90 90]
    fd_sca = fd;                               % Doppler shift
    
    Scample_rate_1 = 32500000;        % scamping rate_1: 32.5MHZ
    Scample_rate_2 = 16250000;        % scamping rate_2: 16.25MHZ
    Scample_rate_3 = 8125000;         % scamping rate_3: 8.125MHZ
    Scample_rate_4 = 4062500;         % scamping rate_4: 4.0625MHZ
    Scample_rate_5 = 2.03125;         % scamping rate_5: 2.03125MHZ
    
                                      % time delay for scatter wave convert time delay point
    delaypoint_sca =  floor(time_delay_sit.*Scample_rate_1);    
    
    A_c = 300000000;                  % speed of light
    lambda = A_c/A_fc;    
    load('location.mat');             % load阵元xyz坐标值
    location = location/200*2.5;      %% 阵元间距与波长比直接影响到谱峰搜索真实DOA的结果，取值不当真实DOA取值结果不如虚假峰值
    order  = 149;   
    Hd = filter_30M_70_6_5m;               % Fir filter; order=16 scampleing rate=32.5M band pass=5MHZ 
    
    N  = 32768*2;                             % the length of generated IQ squence
    N2 = 32768;                               % the length of IQ squence Decimated by 2 :2倍抽取
    Power_noise =0.5*1.38*(1e-23)*(10e6)*1*290;
    Pn = 10*log10(Power_noise)+30;
    Noise_sig = wgn(8,N2,Pn,'dBm','complex'); % gengerate Noise 
    SNR_dir = 10*log10(P_dir/(Power_noise));
    SNR_sca = 10*log10(P_sca/(Power_noise)); 

%%                                      情景2 场景布置图
   figure
    axis([-700 700 -250 250]); hold on;
    scatter(A_Txyz(1,:),A_Txyz(2,:),'MarkerEdgeColor','b', 'LineWidth',1.5);hold on; 
    text(A_Txyz(1,:),A_Txyz(2,:)-10,'Target');
    scatter(A_Sxyz(1,1:A_sca_num), A_Sxyz(2,1:A_sca_num),'d','MarkerEdgeColor','k', 'LineWidth',1.5);hold on; 
    text(A_Sxyz(1,1),A_Sxyz(2,1)-10,'S1'); text(A_Sxyz(1,2),A_Sxyz(2,2)-10,'S2'); text(A_Sxyz(1,3),A_Sxyz(2,3)-10,'S3');
    legend('目标位置','散射体位置');
    scatter(A_Rxyz(1),A_Rxyz(2),'*','MarkerEdgeColor','k','MarkerFaceColor',[0,0.5,0.5]); hold on;
    text(A_Rxyz(1),A_Rxyz(2)-10,'R');
    xlabel('m/米');
    ylabel('m/米');
    title('基于升空散射体计算服务器情景2仿真场景')
    grid minor;
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
        dir_signal_N2 = dir_signal_N(1:decimation:N2);                           % IQ squence decimated by 2  length:32768    
        
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
            attenuation_sca(i_d,i_s)= 10^((SNR_sca(i_d,i_s)-SNR_dir(i_d))/20)*exp(-1i*(A_fc+fd_sca(i_d,i_s))*delaypoint_sca(i_d,i_s)/Scample_rate_1);          % Computer factor: 10*log10(A^2)=SNR_sca(m)-SNR_dir(m)
            tao_sca = zeros(1,8);                                               % 散射波阵元和阵心的信号间隔时延
            sca_array_temple = zeros(8,N2); 
         
            sca_signal_N(1+delaypoint_sca(i_d,i_s):end) = dir_signal_N(1:N-delaypoint_sca(i_d,i_s));  % gengerate scatter wave through direct wave 
            sca_signal_N2=sca_signal_N(1:decimation:N2);                                     % decimate by 2 for scatter wave   两倍抽取 时延点为原来的两倍
            sca_record_array(i_s,:)=attenuation_sca(i_d,i_s)*sca_signal_N2(1,:).*exp(1i*fd_sca(i_d,i_s)*([0:N2-1]/Scample_rate_1*2));      
            count = 0;
            for i = 1:2:16
               tao_sca(count+1) = (location(i,1)*cosd(azimuth_sca(i_s))*sind(pangle_sca(i_s))+location(i,2)*sind(azimuth_sca(i_s))*sind(pangle_sca(i_s))+location(i,3)*cosd(pangle_sca(i_s)))/A_c;
               sca_array_temple(count+1,:) = attenuation_sca(i_d,i_s)*sca_signal_N2(1,:).*exp(1i*fd_sca(i_d,i_s)*([0:N2-1]/Scample_rate_1*2)).*exp(+1i*2*pi*A_c*tao_sca(count+1)/lambda);
               count = count+1;
            end
            sca_array_sig = sca_array_sig+sca_array_temple;
        end
        Array_sig=Array_sig+dir_array_sig+sca_array_sig;   
    end
        Array_sig = Array_sig+Noise_sig;    

 %% DOA estimate 
      Rxx_doa = (1/N2)*(Array_sig*Array_sig');    % generate Covariance matrix about Signal
      [Rdoa_fv,Rdoa_fval]=eig(Rxx_doa);           % 特征分解
      Rdoa_fval2=diag(Rdoa_fval)';                % 特征值向量化
      Rdoa_fval2=abs(Rdoa_fval2);                 
      [Rdoa_fval3,Rdoa_I]=sort(Rdoa_fval2,'descend'); % 特征值从大到小排序
      Rdoa_fv2=Rdoa_fv(1:end,Rdoa_I);                 % 特征向量排序     
      
      Rdoa_fval3_2 = Rdoa_fval3./(min(Rdoa_fval3));  % 特征向量归一化除以最小值
      sig_div1 = zeros(1,7);                         
      for i=1:7
          sig_div1(i)= Rdoa_fval3_2(i)/Rdoa_fval3_2(i+1);
      end
      for i=1:7
          if sig_div1(i)<3
              sig_num=i-1;                           % 指标4估计的信号源数: 直达波数、散射波数(升空散射体可被测向)
              break;
          end
      end
  %%
   Rdoa_Us = Rdoa_fv2(:,1:sig_num);              %% Us信号子空间矩阵  包含直达波和散射波
   Rdoa_Un = Rdoa_fv2(:,sig_num+1:end);          %% Un噪声子空间矩阵  只包含噪声
   
   doa_azimuth =1:0.5:360;
   doa_pangle = 1:1:89;
   doa_array = zeros(size(doa_azimuth,2),size(doa_pangle,2));

   count = 0; am=zeros(8,1);am_dir=zeros(8,1);tao_expectdir=zeros(8,1);
   for i=1:size(doa_azimuth,2)
       for j=1:size(doa_pangle,2)
           for k = 1:2:16
               tao_dir(i) = (location(k,1)*cosd(doa_azimuth(i))*sind(doa_pangle(j))+location(k,2)*sind(doa_azimuth(i))*sind(doa_pangle(j))+location(k,3)*cosd(doa_pangle(j)))/A_c;
               am(count+1)=exp(+1i*2*pi*A_c*tao_dir(i)/lambda);
               count = count+1;
           end
           count =0;
           doa_array(i,j)=1/(am'* (Rdoa_Un* Rdoa_Un')*am);
       end
   end
 figure
[Y,X] =meshgrid(1:89,1:0.5:360);
surf(X,Y,20*log10(abs(doa_array)),'Edgecolor','none')
xlabel('azimuth / °')
ylabel('pitchangle / °');
zlabel('Amplitude');

                                          % Doa Spectral Peak Search
cu_scale = 3;  cu_scale2= 2;
doa_roughsearch = cell(1,45*240); count =0;
doasearch_sumarea = zeros(1,45*240);

 for i=1:239                               %% 区域坐标存储 区域峰值存储
     for j=1:44
        doa_roughsearch{count+1} = [doa_azimuth(cu_scale*(i-1)+1:cu_scale*(i-1)+cu_scale),doa_pangle(cu_scale2*(j-1)+1:cu_scale2*(j-1)+cu_scale2)];
        doasearch_sumarea(count+1) = sum(sum(doa_array(cu_scale*(i-1)+1:cu_scale*(i-1)+cu_scale,cu_scale2*(j-1)+1:cu_scale2*(j-1)+cu_scale2),2));
        count =count+1;
     end
 end
[doasearch_sumarea2,doa_areaIndex] = sort(doasearch_sumarea,'descend');
buchang =10;
doa_detailsearch =cell(1,A_dir_num+A_sca_num+buchang);
for i=1:A_dir_num+A_sca_num+buchang
    doa_detailsearch{i}=doa_roughsearch{doa_areaIndex(i)};
end

 count=0;     %% detail search doa martix

doa_peak = zeros(1,(A_dir_num+A_sca_num)+buchang);                    %% End find peak
doa_endfind =cell(1,(A_dir_num+A_sca_num)+buchang);                   %% pre End find doa 
doa_endazmiuth =zeros(1,(A_dir_num+A_sca_num)+buchang);
doa_endfind2 =cell(1,(A_dir_num+A_sca_num));                          %% End find doa
for i=1:A_dir_num+A_sca_num+buchang
    doa_peak_temp = zeros(1,200);
    doa_detail = cell(1,200); 
   for j=doa_detailsearch{i}(1):0.1:doa_detailsearch{i}(cu_scale)
       for k=doa_detailsearch{i}(cu_scale+1):0.1:doa_detailsearch{i}(cu_scale+cu_scale2)
           doa_detail{count+1} =[j,k];
           count2=0;
           for e = 1:2:16
               tao_dir_temple = (location(e,1)*cosd(j)*sind(k)+location(e,2)*sind(j)*sind(k)+location(e,3)*cosd(k))/A_c;
               am(count2+1)=exp(+1i*2*pi*A_c*tao_dir_temple/lambda);
               count2 = count2+1;
           end
           doa_peak_temp(count+1)=1/(am'* (Rdoa_Un* Rdoa_Un')*am);
           count =count+1;
       end
   end
   [doa_peak_value,peak_index]=max(doa_peak_temp);
    doa_peak(i)= doa_peak_value;
    doa_endfind{i} = doa_detail{peak_index};
    doa_endazmiuth(i)=doa_detail{peak_index}(1);         %% 存取方向角
    count =0;
end
   for i=1:size(doa_peak,2)
       temp=doa_endazmiuth(i);
       value = doa_peak(i);
       temp_md =abs(temp-doa_endazmiuth);
       for j=1:size(doa_peak,2)
           if(temp_md(j)<=4 && temp_md(j)>0 )
               if doa_peak(j)<=value
                   doa_peak(j)=0;
                   doa_endazmiuth(j)=0;
               else
                   doa_peak(i)=0;
                   doa_endazmiuth(i)=0;
                   break;
               end
           end
       end
   end
    [~,endIndex]=sort(doa_peak,'descend');
    for i=1:A_dir_num+A_sca_num
        doa_endfind2{i}=doa_endfind{endIndex(i)};
    end

 %% 参考通道波束形成 与 扇区波束形成 
   sig_resource = 3;                         %% Recording the number of current signal by DOA_Estimate_argrithm 
   sig_estimate_resource =1;                 %% signal resource by DOA_estimate_algorithm
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
   if (sig_resource > sig_estimate_resource)                                  % Computer Service situtation 2
       for i = 1:2:16                                                         % default 1 direct wave
           tao_dir(i) = (location(i,1)*cosd(azimuth_dir(1))*sind(pangle_dir(1))+location(i,2)*sind(azimuth_dir(1))*sind(pangle_dir(1))+location(i,3)*cosd(pangle_dir(1)))/A_c;
           vectorA_ref(count+1)=exp(+1i*2*pi*A_c*tao_dir(i)/lambda);
           count = count+1;
       end
       Wref = (Rx_Msca\vectorA_ref(:,1))/(vectorA_ref(:,1)'/Rx_Msca*vectorA_ref(:,1));
       y_ref = Wref'*Array_sig;
   elseif(sig_resource == sig_estimate_resource)                              %% Computer Service situtation 3                                                                    %% Computer Service situtation 3
       if A_dir_num>1
       [~,index_sca]=max(max(SNR_sca));
      else
       [~,index_sca]=max(SNR_sca);   
      end
       tao_sca  = zeros(1,8); count=0;
       for i=1:2:16
           tao_sca(count+1) = (location(i,1)*cosd(azimuth_sca(index_sca))*sind(pangle_sca(index_sca))+location(i,2)*sind(azimuth_sca(index_sca))*sind(pangle_sca(index_sca))+location(i,3)*cosd(pangle_sca(index_sca)))/A_c;
           count = count+1;
       end
       vectorA_ref=exp(+1i*2*pi*A_c*tao_sca'/lambda);
       Wref = (Rx_Msca\vectorA_ref(:,1))/(vectorA_ref(:,1)'/Rx_Msca*vectorA_ref(:,1));
       y_ref = Wref'*Array_sig;
   end
   
 %%% Generate 12 Sector beam  BLOCK 
   tao_sca  = zeros(1,8);
   Sector_array=[0 30 60 90 120 150 180 210 240 270 300 330;
                10 10 10 10  10  10  10  10  10  10  10  10];
   VectorA_sector = zeros(8,12);                                    %% Array vector for every sector(扇区)
   VectorA_sca = zeros(8,A_sca_num); 
   Wsector_array_mvdr = zeros(8,12);                                %% computer w based MVDR
   Wsector_array_block = zeros(8,12);                               %% computer w based Block array
   Wscatter_array = zeros(8,A_sca_num);
   z_sector = zeros(12,N2);
   z_sector_block = zeros(12,N2);
   for i_s =1:12
       count =0;
       for i=1:2:16
           tao_sca(count+1) = (location(i,1)*cosd(Sector_array(1,i_s))*sind(Sector_array(2,i_s))+location(i,2)*sind(Sector_array(1,i_s))*sind(Sector_array(2,i_s))+location(i,3)*cosd(Sector_array(2,i_s)))/A_c;
           count = count+1;
       end
       VectorA_sector(:,i_s) = exp(+1i*2*pi*A_c*tao_sca'./lambda);
       Wsector_array_mvdr(:,i_s) = (Rx_Msca\VectorA_sector(:,i_s))/(VectorA_sector(:,i_s)'/Rx_Msca*VectorA_sector(:,i_s));
       Wsector_array_block(:,i_s) = Rx_Block*VectorA_sector(:,i_s);
   end
   
   for i_s =1:A_sca_num
       count =0;
       for i=1:2:16
         tao_sca(count+1) = (location(i,1)*cosd(azimuth_sca(i_s))*sind(pangle_sca(i_s))+location(i,2)*sind(azimuth_sca(i_s))*sind(pangle_sca(i_s))+location(i,3)*cosd(pangle_sca(i_s)))/A_c;
         count = count+1;
       end
       VectorA_sca(:,i_s) = exp(+1i*2*pi*A_c*tao_sca'./lambda);
   end
   Wscatter_array(:,1) = (Rx_Msca\VectorA_sca(:,1))/(VectorA_sca(:,1)'/Rx_Msca*VectorA_sca(:,1));
   Wscatter_array(:,2) = (Rx_Msca\VectorA_sca(:,2))/(VectorA_sca(:,2)'/Rx_Msca*VectorA_sca(:,2));
   z1 = Wscatter_array(:,1)'*Array_sig;
   z1_2=Wscatter_array(:,2)'*Array_sig;
   
   for i=1:12
       z_sector(i,:)=Wsector_array_mvdr(:,i)'*Array_sig;                    %% Sector beamforming based MVDR
       z_sector_block(i,:)=Wsector_array_block(:,i)'*Array_sig;             %% Sector beamforming based block
   end         

%% Generate RD spectrum
   scampling =Scample_rate_1;
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
   window_h =(0.35875*ones(1,N2)-0.48829*cos((2*pi*[1:N2])/N2)+0.1428*cos((4*pi*[1:N2])/N2)-0.01168*cos((6*pi*[1:N2])/N2));
   RD_cell=cell(1,12);
   k=1;
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
[Y,X] =meshgrid(-k:k,-l:l);
Rd_plot = RD_cell{1,1};
surf(X,Y,abs(Rd_plot),'Edgecolor','none')
xlabel('time delay / point')
ylabel('Dopple frequence / point');
zlabel('Amplitude/dB');
%% CFAR检验
if scampling ==Scample_rate_1
    N_pf=200;  L_k = [163,1];
elseif scampling ==Scample_rate_2
    N_pf=150;  L_k = [82,2];
elseif scampling ==Scample_rate_3
    N_pf=70;
elseif scampling ==Scample_rate_4
    N_pf=30;
elseif scampling ==Scample_rate_5
    %l = 1;  k=113;
end
Rd_array=RD_cell{1,1};
feng = zeros(1,L_k(2));
for j=1:2*L_k(2)+1
    feng(j)=sum(abs(Rd_array(:,j)));   % 选取峰值所在列
end
[~,cfar_index]=max(feng);
Rd_sector = abs(RD_cell{1,10});                                          % 取出第一扇区的数据
flag_noise = NoisType(Rd_sector(:,cfar_index)',0.0001);     
disp(['噪声类型：',sprintf('%d',flag_noise)])
Rd_date = abs(Rd_array(:,cfar_index));
Pf=10e-4; l =2*163+1;  N=20;  k=floor(N*3/4);  Product_unit = 2;        % Rd_数据必须是行序列
[Order_L,T_Sn,xk]=Qin_CFAR(Pf,Rd_date',l,flag_noise,N_pf,N,k,Product_unit); % l+delay_point 应该是CFAR输出的对应delya_point点
disp(['CFAR门限值：',sprintf('T:%0.2f, Sn:%0.2f',T_Sn(1),T_Sn(2))])
%%  解析帧文件生成
frame_decoder = load('E:\混合数据帧2解析\Mixframe2_2020_11_23_15_49_42_493_decode.mat');
frame_decoder = frame_decoder.Struct_data;
frame_decoder{20} = size(A_Sxyz,2);  % Num_Q
Num_T = 0;
frame_decoder{56} = Num_T;           % Num_T
Num_DOA =4;                                                                   % 直达波数量 + 散射波数量
frame_decoder{41} = Num_DOA;         % Num_DOA
% 2目标
% DOA = [azimuth_dir(1) pitch_dir_sit(1)                         
%        azimuth_dir(2) pitch_dir_sit(2)                                         
%        azimuth_sca(3) pitch_sca_sit(3) 
%        azimuth_sca(1) pitch_sca_sit(1) 
%        azimuth_sca(2) pitch_sca_sit(2)   
%         0     0
%         0     0 ];
% 1目标
DOA = [azimuth_dir(1) pitch_dir_sit(1)                                                                 
       azimuth_sca(3) pitch_sca_sit(3) 
       azimuth_sca(1) pitch_sca_sit(1) 
       azimuth_sca(2) pitch_sca_sit(2)   
        0     0
        0     0 
        0     0];
% 2目标
% T_info = [1,4                                                   
%           2,5
%           3,3
%           0, 0];
% 1目标
T_info = [0,0                                                   
          0,0
          0,0
          0, 0];
frame_decoder{1,65} = 1;                % 参考散射体的DOA序号
inital = 1;
for i=1:7
    frame_decoder{1,41+inital} =  DOA(i,1);                                  % 方向角0-360
    frame_decoder{1,41+inital+1}= DOA(i,2);                                  % 俯仰角0-90
    inital = inital + 2;
end
Sxyz = [A_Sxyz(1,1),A_Sxyz(2,1),A_Sxyz(3,1),0,azimuth_sca(1);
       A_Sxyz(1,2),A_Sxyz(2,2),A_Sxyz(3,2),0,azimuth_sca(2);
      A_Sxyz(1,3) ,A_Sxyz(2,3),A_Sxyz(3,3),0,azimuth_sca(3);
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
frame_order = 5;
str = sprintf('E:\\混合数据帧2解析\\Mixframe2_2020-11-Qin2-目标个数%d-序号%d_decode',A_dir_num,frame_order);
save(str,'Struct_data')   
%%
%  %%  仿真图
%      xaxi = linspace(-N2+1,N2-1,2*N2-1);   
%      %%% 绘制阵列信号 对应的功率谱图
%      figure
%      w1 = window('hamming',N2);
%      fs1 = 16250000;
%      f1 = -fs1/2:fs1/N2:fs1/2-1;           %%频谱分辨率
%      S_fft=Array_sig(1,1:N2).* w1';
%      plot(f1,20*log10(abs(fftshift(fft(S_fft)))));
%      xlabel('HZ');
%      ylabel('Amplitude/dB');
%      title ('8通道阵列接收设备 载波频率:1.5GHZ 带宽:10MHZ ')
%      grid on
%      %%% 绘制阵列信号 自相关图 验证直达波和散射波是否存在 散射波的位置是否正确
%      figure
%      cor_S_1 =20*log10(abs(xcorr(Array_sig(1,1:N2))));
%      plot(xaxi,cor_S_1)
%      title('测试： 阵列信号自相关')
% %      set(gca, 'GridAlpha', 1,'GridColor',[0 0 0],'MinorGridColor',[0 0 0],'XMinorTick','on','YMinorTick','on');  % 设置透明度
%      xlabel('time delay / point')
%      ylabel('Amplitude/dB');
%      grid on
% 
% %%
%    corrcoef(dir_signal_N2(1,:),    y_ref)               %%计算协方差系数
%    corrcoef(sca_record_array(1,:), y_ref)               %%计算协方差系数
% %  corrcoef(sca_record_array(2,:), z1_2)                %%计算协方差系数
%    figure
%    subplot(2,1,1)
%      w1 = window('hamming',N2);
%      fs1 = 16250000;
%      f1 = -fs1/2:fs1/N2:fs1/2-1;           %%频谱分辨率
%      S_fft=dir_signal_N2(1,1:N2).* w1';
%      plot(f1,20*log10(abs(fftshift(fft(S_fft)))));
%      xlabel('HZ');
%      ylabel('Amplitude/dB');
%      title ('8通道阵列接收设备 直达波信号频谱图 带宽:10MHZ ')
%      grid on
%    subplot(2,1,2)
%        w1 = window('hamming',N2);
%        fs1 = 16250000;
%        f1 = -fs1/2:fs1/N2:fs1/2-1;           %%频谱分辨率
%        S_fft=y_ref(1,1:N2).* w1';
%        plot(f1,20*log10(abs(fftshift(fft(S_fft)))));
%        xlabel('HZ');
%        ylabel('Amplitude/dB');
%        title ('8通道阵列接收设备 参考通道信号频谱图 带宽:10MHZ ')
%        grid on
%    figure
%     subplot(3,1,1)
%      cor_S_1 =20*log10(abs(xcorr(Array_sig(1,1:N2),y_ref)));
%      plot(xaxi,cor_S_1)
%      title('测试： 阵列信号与参考通道信号互相关')
% %      set(gca, 'GridAlpha', 1,'GridColor',[0 0 0],'MinorGridColor',[0 0 0],'XMinorTick','on','YMinorTick','on');  % 设置透明度
%      xlabel('time delay / point')
%      ylabel('Amplitude/dB');
%      grid on
%     subplot(3,1,2)
%      cor_S_2 =20*log10(abs(xcorr(y_ref,z_sector_block(1,:))));
%      plot(xaxi,cor_S_2)
%      title('测试：参考通道信号与波束形成：扇区2信号互相关')
% %      set(gca, 'GridAlpha', 1,'GridColor',[0 0 0],'MinorGridColor',[0 0 0],'XMinorTick','on','YMinorTick','on');  % 设置透明度
%      xlabel('time delay / point')
%      ylabel('Amplitude/dB');
%      grid on
%     subplot(3,1,3)
%      cor_S_3 =20*log10(abs(xcorr(y_ref,z_sector_block(7,:))));
%      plot(xaxi,cor_S_3)
%      title('测试： 参考通道信号与波束形成：扇区信号7互相关')
% %      set(gca, 'GridAlpha', 1,'GridColor',[0 0 0],'MinorGridColor',[0 0 0],'XMinorTick','on','YMinorTick','on');  % 设置透明度
%      xlabel('time delay / point')
%      ylabel('Amplitude/dB');
%      grid on
%   figure
%      cor_S_3 =20*log10(abs(xcorr(y_ref,z_sector_block(10,:))));
%      plot(xaxi,cor_S_3)
%      title('测试： 参考通道信号与波束形成：扇区信号10互相关')
% %      set(gca, 'GridAlpha', 1,'GridColor',[0 0 0],'MinorGridColor',[0 0 0],'XMinorTick','on','YMinorTick','on');  % 设置透明度
%      xlabel('time delay / point')
%      ylabel('Amplitude/dB');
%      grid on
% %%
%    similar_sectorwave = zeros(1,12);
%    for i=1:A_sca_num
%        for j=1:12
%            cor = corrcoef(sca_record_array(i,:), z_sector_block(j,:));                    %%计算协方差系数
%            similar_sectorwave(i,j)=abs(cor(1,2));
%        end
%    end
%    figure
%    subplot(2,1,1)
%      w1 = window('hamming',N2);
%      fs1 = 16250000;
%      f1 = -fs1/2:fs1/N2:fs1/2-1;           %%频谱分辨率
%      S_fft=sca_record_array(1,1:N2).* w1';
%      plot(f1,20*log10(abs(fftshift(fft(S_fft)))));
%      xlabel('HZ');
%      ylabel('Amplitude/dB');
%      title ('8通道阵列接收设备 方位角13°散射波信号频谱图 带宽:10MHZ ')
%      grid on
%    subplot(2,1,2)
%        w1 = window('hamming',N2);
%        fs1 = 16250000;
%        f1 = -fs1/2:fs1/N2:fs1/2-1;         %%频谱分辨率
%        S_fft=z_sector_block(5,1:N2).* w1';
%        plot(f1,20*log10(abs(fftshift(fft(S_fft)))));
%        xlabel('HZ');
%        ylabel('Amplitude/dB');
%        title ('8通道阵列接收设备 0°扇区1波束形成信号频谱图 带宽:10MHZ ')
%        grid on  
%    figure
%    subplot(2,1,1)
%      w1 = window('hamming',N2);
%      fs1 = 16250000;
%      f1 = -fs1/2:fs1/N2:fs1/2-1;           %%频谱分辨率
%      S_fft=sca_record_array(2,1:N2).* w1';
%      plot(f1,20*log10(abs(fftshift(fft(S_fft)))));
%      xlabel('HZ');
%      ylabel('Amplitude/dB');
%      title ('8通道阵列接收设备  方位角170°散射波信号频谱图 带宽:10MHZ ')
%      grid on
%    subplot(2,1,2)
%        w1 = window('hamming',N2);
%        fs1 = 16250000;
%        f1 = -fs1/2:fs1/N2:fs1/2-1;           %%频谱分辨率
%        S_fft=z_sector_block(7,1:N2).* w1';
%        plot(f1,20*log10(abs(fftshift(fft(S_fft)))));
%        xlabel('HZ');
%        ylabel('Amplitude/dB');
%        title ('8通道阵列接收设备 扇区7波束形成信号频谱图 带宽:10MHZ ')
%        grid on  
%    figure
%    subplot(2,1,1)
%      w1 = window('hamming',N2);
%      fs1 = 16250000;
%      f1 = -fs1/2:fs1/N2:fs1/2-1;           %%频谱分辨率
%      S_fft=sca_record_array(3,1:N2).* w1';
%      plot(f1,20*log10(abs(fftshift(fft(S_fft)))));
%      xlabel('HZ');
%      ylabel('Amplitude/dB');
%      title ('8通道阵列接收设备 方位角270°直达波信号频谱图 带宽:10MHZ ')
%      grid on
%    subplot(2,1,2)
%        w1 = window('hamming',N2);
%        fs1 = 16250000;
%        f1 = -fs1/2:fs1/N2:fs1/2-1;           %%频谱分辨率
%        S_fft=z_sector_block(10,1:N2).* w1';
%        plot(f1,20*log10(abs(fftshift(fft(S_fft)))));
%        xlabel('HZ');
%        ylabel('Amplitude/dB');
%        title ('8通道阵列接收设备 扇区10参考通道信号频谱图 带宽:10MHZ ')
%        grid on  
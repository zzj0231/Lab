A_Txyz = [400,120,5]';  nf =0.78; sig_num=3; 
A_Sxyz = [600,65,50;
        -700,70,60]';
sim_value=zeros(2,100);
rmse_az=zeros(1,16);
rmse_pt=zeros(1,16);
index=1;
all_nf=[1.95,1.65,1.25,0.95,0.78,0.66,0.52,0.41,0.33,0.25,0.2,0.16,0.12,0.1,0.08,0.016];
for ii=1:16
    nf=all_nf(ii);
    for i=1:100
        [Array_sig,doaStandDir,doaStandSca,SNR_dir,SNR_sca,delaypoint_sca]=produceArraySig(A_Txyz,A_Sxyz,nf);
        [doa_roughsearch2,doa_EndFind] = estiamteDoa(Array_sig,sig_num);
        index_doa = find_doa(doaStandSca(:,1),doa_EndFind);
        if (index_doa~=0)
            sim_value(:,i)=doa_EndFind{index_doa}';
        end
    end
%     static(ii)=size(sim_value,2);
    sim_value(:,all(sim_value==0,1))=[];
    rmse_az(index) = RMSE(doaStandSca(1,1),sim_value(1,:));
    rmse_pt(index) = RMSE(doaStandSca(2,1),sim_value(2,:));
    index=index+1;
end





function [Array_sig,doaStandDir,doaStandSca,SNR_dir,SNR_sca,delaypoint_sca] = produceArraySig(A_Txyz,A_Sxyz,nf)
location=load('location.mat');             
location = location.location/200*2.5;      %% 阵元间距与波长比直接影响到谱峰搜索真实DOA的结果，取值不当真实DOA取值结果不如虚假峰值
% Gengerate IQ data and Process IQ data FOR COMPUTER SECOND SEVERCE
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
    fd_sca = fd;                            
    
    Scample_rate_1 = 32500000;        % scamping rate_1: 32.5MHZ
    
                                      % time delay for scatter wave convert time delay point
    delaypoint_sca =  floor(time_delay_sit.*Scample_rate_1);    
    
    A_c = 300000000;                  % speed of light
    lambda = A_c/A_fc;    
    
    order  = 149;   
    Hd = filter_65M_70_6_5m;               % Fir filter; order=16 scampleing rate=32.5M band pass=5MHZ 
    
    N  = 32768*3;                             % the length of generated IQ squence
    N2 = 32768;                               % the length of IQ squence Decimated by 2 :2倍抽取
    
    Power_noise =nf*1.38*(1e-23)*(30e6)*1*290;
    Pn = 10*log10(Power_noise)+30;
    Noise_sig = wgn(8,N2,Pn,'dBm','complex'); % gengerate Noise 
    SNR_dir = 10*log10(P_dir/(Power_noise));
    SNR_sca = 10*log10(P_sca/(Power_noise)); 
    
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
end

function [doa_roughsearch2,doa_EndFind] = estiamteDoa(Array_sig,sig_num)
  location=load('location.mat');
  location = location.location/200*2.5;
  A_c = 300000000;                  % speed of light
  A_fc = 1500000000;                  % carrier frequence:1.5G
  lambda = A_c/A_fc;
  N2=32678;
  
  Rxx_doa = (1/N2)*(Array_sig*Array_sig');   
  [Rdoa_fv,Rdoa_fval]=eig(Rxx_doa);          
  Rdoa_fval2=diag(Rdoa_fval)';                % 特征值向量化
  Rdoa_fval2=abs(Rdoa_fval2);
  [~,Rdoa_I]=sort(Rdoa_fval2,'descend');      % 特征值从大到小排序
  Rdoa_fv2=Rdoa_fv(1:end,Rdoa_I);            
  Rdoa_Un = Rdoa_fv2(:,sig_num+1:end);       
   
   tao_dir=zeros(1,8); 
% 粗搜索范围
plotazimuth =0:2:360;
plotpangle = 0:1:89;
plotarray = zeros(size(plotazimuth,2),size(plotpangle,2));

count = 0; am=zeros(8,1);
for i=1:size(plotazimuth,2)
   for j=1:size(plotpangle,2)
       for k = 1:2:16
           tao_dir(i) = (location(k,1)*cosd(plotazimuth(i))*sind(plotpangle(j))+location(k,2)*sind(plotazimuth(i))*sind(plotpangle(j))+location(k,3)*cosd(plotpangle(j)))/A_c;
           am(count+1)=exp(+1i*2*pi*A_c*tao_dir(i)/lambda);
           count = count+1;
       end
       count =0;
       plotarray(i,j)=1/(am'* (Rdoa_Un* Rdoa_Un')*am);
   end
end
[music,fangwei,fuyang]=peak_search(plotarray,0,0,2,1);  % 修改处 如果粗搜索的步长发生改变
doa_roughsearch2 = cell(1,size(fangwei,2)); 
[~,index]=sort(music,'descend');
for i=1:size(fangwei,2)
    doa_roughsearch2{i}=[fangwei(index(i)),fuyang(index(i))];
end
doa_detailsearch2 = cell(1,size(fangwei,2)); 
for i=1:size(doa_roughsearch2,2)
    az = doa_roughsearch2{i}(1);
    pt = doa_roughsearch2{i}(2);
    az_min = az-3; az_max=az+3; pt_min = pt-1; pt_max=pt+1;
    if az<3
        az_min=1;
    elseif az>357
         az_max=360;   
    end
    if pt<1
        pt_min =1;
    elseif pt>89
        pt_max = 90;
    end
    doa_detailsearch2{i}=[az_min,az_max,pt_min,pt_max];
end
doa_EndFind = detailsearch(doa_detailsearch2,0.1,Rdoa_Un,location,1.5e9);
end

function rmse = RMSE(ref_value,sim_value)
   N = size(sim_value,2);
   rmse = sqrt(sum(sum((sim_value-ref_value).^2))*1/N);
end

function index = find_doa(ref_doa,sim_value)
   N = size(sim_value,2);
   index=0;
   doa=zeros(2,N);
   for i=1:N
      doa(1,i)=sim_value{i}(1);
      doa(2,i)=sim_value{i}(2);
   end
  [error,index_order] = min(sqrt(sum((ref_doa-doa).^2)));
   if error<3
       index=index_order;
   end
end
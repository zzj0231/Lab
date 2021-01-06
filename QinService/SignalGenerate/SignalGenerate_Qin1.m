% Gengerate IQ data and Process IQ data FOR COMPUTER SECOND SEVERCE
    A_Rxyz = [5,-10,5]';       % location of Receiver column 
    A_Txyz = [15,180,5]';    % location of Target   column
    A_Sxyz = [-25,200,15;
               60,200,35]'; 
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
    Scample_rate_5 = 2.03125;          % scamping rate_5: 2.03125MHZ
    
                                      % time delay for scatter wave convert time delay point
    delaypoint_sca =  floor(time_delay_sit.*Scample_rate_1);    
    
    A_c = 300000000;                  % speed of light
    lambda = A_c/A_fc;    
    load('location.mat');             % load阵元xyz坐标值
    location = location/200*2.5;      %% 阵元间距与波长比直接影响到谱峰搜索真实DOA的结果，取值不当真实DOA取值结果不如虚假峰值
    order  = 149;   
    Hd = filter_65M_70_6_5m;               % Fir filter; order=16 scampleing rate=32.5M band pass=5MHZ 
    
    N  = 32768*3;                             % the length of generated IQ squence
    N2 = 32768;                               % the length of IQ squence Decimated by 2 :2倍抽取
    Power_noise =0.17*1.38*(1e-23)*(30e6)*1*290;
    Pn = 10*log10(Power_noise)+30;
    Noise_sig = wgn(8,N2,Pn,'dBm','complex'); % gengerate Noise 
    SNR_dir = 10*log10(P_dir/(Power_noise));
    SNR_sca = 10*log10(P_sca/(Power_noise)); 

%%                                    情景1 单站轨迹图
A_RxyzTrace=[-60,5,5;
              5,-10,5
              75,10,5]';
az_change = zeros(1,size(A_RxyzTrace,2)-1);
for i=1:size(az_change,2)
    R1=A_RxyzTrace(:,i);
    R2=A_RxyzTrace(:,i+1);
    az_change(i)=atan((R2(2)-R1(2))/(R2(1)-R1(1)))*180/pi;
end
   figure
    axis([-80 100 -30 350]); hold on;
    scatter(A_Txyz(1,:),A_Txyz(2,:),'MarkerEdgeColor','b', 'LineWidth',1.5);hold on; 
    text(A_Txyz(1,:)+5,A_Txyz(2,:)-5,'Target');
    scatter(A_Sxyz(1,1:A_sca_num), A_Sxyz(2,1:A_sca_num),'v','MarkerEdgeColor','k', 'LineWidth',1.5);hold on; 
    text(A_Sxyz(1,1),A_Sxyz(2,1)-13,'S1'); text(A_Sxyz(1,2),A_Sxyz(2,2)-13,'S2'); %text(A_Sxyz(1,3),A_Sxyz(2,3)-13,'S3');
    
    scatter(A_RxyzTrace(1,:),A_RxyzTrace(2,:),'^','MarkerEdgeColor','k','MarkerFaceColor',[0,0.5,0.5]); hold on;
    text(A_Rxyz(1),A_Rxyz(2)-15,'R');
    legend('目标位置','反射物位置','地面单站');
    
    xlabel('m/米');
    ylabel('m/米');
    title('移动单站目标定位-oxy')
    set(gca, 'GridAlpha', 1,'GridColor',[0 0 0],'MinorGridColor',[0 0 0],'XMinorTick','on','YMinorTick','on');  % 设置透明度
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
          sig_div1(i)= Rdoa_fval3_2(i)-Rdoa_fval3_2(i+1);
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

   count = 0; am=zeros(8,1);
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
[Y,X] =meshgrid(1:89,1:0.5:360);
surf(X,Y,(abs(doa_array)/10),'Edgecolor','none')
xlabel('方位角 / °')
ylabel('俯仰角 / °');
zlabel('幅度');
title('DOA谱峰搜索')
% 粗搜索范围
plotazimuth =0:2:360;
plotpangle = 0:1:89;
plotarray = zeros(size(plotazimuth,2),size(plotpangle,2));

count = 0; am=zeros(8,1);am_dir=zeros(8,1);tao_expectdir=zeros(8,1);
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
figure
axis([0 360 0 89])
[Y,X] =meshgrid(0:1:89,0:2:360);       % 修改处 如果粗搜索的步长发生改变
set(gca,'XTick',0:20:358);  hold on
set(gca,'YTick',0:10:89);  hold on
surf(X,Y,(20*log(abs(plotarray))/10),'Edgecolor','none')   
xlabel('方位角 / °')
ylabel('俯仰角 / °');
zlabel('幅度');
title('DOA谱峰搜索')
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
                                          % Doa Spectral Peak Search
cu_scale = 3;  cu_scale2= 2;
doa_roughsearch = cell(1,45*240); count =0;
doasearch_sumarea = zeros(1,45*240);
buchang =10;
 for i=1:239                               %% 区域坐标存储 区域峰值存储
     for j=1:44
       doa_roughsearch{count+1} = [doa_azimuth(cu_scale*(i-1)+1:cu_scale*(i-1)+cu_scale),doa_pangle(cu_scale2*(j-1)+1:cu_scale2*(j-1)+cu_scale2)];
       doasearch_sumarea(count+1) = max(max(doa_array(cu_scale*(i-1)+1:cu_scale*(i-1)+cu_scale,cu_scale2*(j-1)+1:cu_scale2*(j-1)+cu_scale2)));  % 区间最大值代表全区间 
       count =count+1;
     end
 end

[doasearch_sumarea2,doa_areaIndex] = sort(doasearch_sumarea,'descend');
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
    for i=1:7
        doa_endfind2{i}=doa_endfind{endIndex(i)};
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
   RD_cell=cell(1,12);
   for i=1:12
      Rd_array = zeros(2*l+1,2*k+1);
      RD_cell{1,i} = Rd_array;
   end
%%  解析帧文件生成
frame_decoder = load('Mixframe2_2020_11_23_15_49_42_493_decode.mat');
frame_decoder = frame_decoder.Struct_data;
frame_decoder{12}=A_Rxyz(1);                        % 修改5 若站位发生改变需要修改坐标
frame_decoder{13}=A_Rxyz(2);
frame_decoder{14}=A_Rxyz(3);
frame_decoder{20} = 0;  % Num_Q
Num_T = 0;
frame_decoder{56} = Num_T;           % Num_T
Num_DOA =sig_num;                                   % 修改1 若sig_num估计不准                        
frame_decoder{41} = Num_DOA;         % Num_DOA      

% 1目标                                          
DOA = zeros(7,2);
for i=1:Num_DOA
    DOA(i,1) = doa_EndFind{i}(1,1);
    DOA(i,2) = doa_EndFind{i}(1,2);
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
% 散射体模板                                          
% Sxyz = [A_Sxyz(1,1),A_Sxyz(2,1),A_Sxyz(3,1),0,azimuth_sca(1);
Sxyz = zeros(4,5)';  
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
frame_order = 3;                                     % 修改2 若保存序号发生改变
stat_order = 2;                                      % 修改3 若站数序号发生改变
str = sprintf('E:\\混合数据帧2解析\\定位情景1帧解析\\Mixframe2_Qin1-站%d-S%d-T%d-序号%d_decode',stat_order,A_sca_num,A_dir_num,frame_order);
save(str,'Struct_data')   
 %%  计算定位误差  
 A_T = [55,220,5]';
 Loc_error = sqrt((loc_result(1,:)-A_T(1,1)).^2+(loc_result(2,:)-A_T(2,1)).^2);
[Loc_error ,index_loc] = sort(Loc_error);
for i=1:size(Loc_error,2)
   if (Loc_error(i)>15)
       Locindex_error15 = index_loc(i-1);
       errorindex_15 = i-1;
       break;
   end
end
%%  绘制定位结果
  figure
  axis([-50 50 160 200]); hold on;
  scatter(A_Txyz(1,:),A_Txyz(2,:),'MarkerEdgeColor','b', 'LineWidth',1.5);hold on; 
  scatter(loc_result(1,:),loc_result(2,:),'c','MarkerEdgeColor','k','MarkerFaceColor','r'); hold on;
  xlabel('m / 米');
  ylabel('m / 米');
  set(gca,'GridColor',[0 0 0],'MinorGridColor',[0 0 0],'XMinorTick','on','YMinorTick','on','XMinorGrid','on');  % 设置透明度
  title('双散射体辅助地面单站定位结果-oxy')
  legend('目标位置','目标定位解算位置');
  error_min = sprintf('e:%0.1fm',Loc_error(1));
  error_15 = sprintf('e:%0.1fm',Loc_error(errorindex_15));
  text(loc_result(1,index_loc(1))-1,loc_result(2,index_loc(1))-1,error_min); hold on
  text(loc_result(1,Locindex_error15)-1,loc_result(2,Locindex_error15)-1,error_15); hold on
  grid on
%% 绘制送入到点迹分选区位置点
color={[1 0 0],[0.8500 0.3250 0.0980],'#0072BD','#77AC30','#4DBEEE'};
for i=1:size(locresult_cell,2)
    locresult = locresult_cell{i};
  axis([-200 100 100 550]); hold on;
  scatter(locresult(1,:),locresult(2,:),'MarkerEdgeColor','w', 'MarkerFaceColor',color{i},'LineWidth',1.5);hold on;
end
xlabel('m / 米');
ylabel('m / 米');
set(gca,'GridColor',[0 0 0],'MinorGridColor',[0 0 0],'XMinorTick','on','YMinorTick','on','XMinorGrid','on');  % 设置透明度
mylegend=['第1个CPI']; mylegend2=['第2个CPI']; mylegend3=['第3个CPI']; mylegend4=['第4个CPI']; mylegend5=['第5个CPI'];
legend(mylegend,mylegend2,mylegend3,mylegend4,mylegend5);
title('双散射体辅助地面单站目标5个CPI时刻定位结果-oxy')
grid on
%%  绘制点迹分选
  pos_real=[15,180,5
            25,190,5
            35,200,5
            45,210,5
            55,220,5]';
  pos_trace1 = postrace_cell{1};
  pos_trace2 = postrace_cell{2};
  pos_trace3 = postrace_cell{3};
  axis([-20 70 160 250]); hold on;
  scatter(pos_real(1,:),pos_real(2,:),'MarkerEdgeColor','b', 'LineWidth',1.5);hold on;
  scatter(pos_trace1(1,:),pos_trace1(2,:),'c','MarkerEdgeColor','k','MarkerFaceColor','r'); hold on;
  scatter(pos_trace2(1,:),pos_trace2(2,:),'c','MarkerEdgeColor','k','MarkerFaceColor','#0072BD'); hold on;
  scatter(pos_trace3(1,:),pos_trace3(2,:),'c','MarkerEdgeColor','k','MarkerFaceColor','#77AC30'); hold on;
  plot(pos_real(1,:),pos_real(2,:),'-', 'LineWidth',1,'Color','k');hold on;
  
  plot(pos_trace1(1,:),pos_trace1(2,:),'-', 'LineWidth',0.3,'Color','k');hold on;
  plot(pos_trace2(1,:),pos_trace2(2,:),'-', 'LineWidth',0.3,'Color','k');hold on;
  plot(pos_trace3(1,:),pos_trace3(2,:),'-', 'LineWidth',0.3,'Color','k');hold on;
  
  text(pos_trace1(1,1)-8,pos_trace1(2,2)-8,'轨迹1'); hold on
  text(pos_trace2(1,1)-7,pos_trace2(2,2)-3,'轨迹2'); hold on
  text(pos_trace3(1,1)-9,pos_trace3(2,2)-9,'轨迹3'); hold on
  mylegend=['真实目标轨迹']; mylegend2=['目标定位结果轨迹1']; mylegend3=['目标定位结果轨迹2']; mylegend4=['目标定位结果轨迹3'];
  legend(mylegend,mylegend2,mylegend3,mylegend4);
  xlabel('m / 米');
  ylabel('m / 米');
  set(gca,'GridColor',[0 0 0],'MinorGridColor',[0 0 0],'XMinorTick','on','YMinorTick','on','XMinorGrid','on');  % 设置透明度
  title('双散射体辅助地面单站目标跟踪滤波结果-oxy')
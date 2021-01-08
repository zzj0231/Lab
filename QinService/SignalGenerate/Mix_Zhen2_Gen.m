
%% 设置路径、load特定的Mixframe_decoder、调用函数生成Mixframe2的bin文件
Num_file = 8;
for i=0:Num_file-1
str_recoder = sprintf('D:\\ZZJ_项目\\升空散射体-计算服务器\\混合数据帧2\\MixZhen2_20210107Qin3_%d.bin',i);
   str_load = sprintf('E:\\混合数据帧2解析\\Mixframe2_20210107Qin3-S3-T1-序号%d_decode.mat',i); 
Mixframe_decoder = load(str_load).Struct_data;
Mixframe_decoder{66} = zeros(8,8);
generate_Mixframe(Mixframe_decoder,str_recoder);
end
%%
fileId =fopen(pathQin3_3,'rb');
mixframe2_bin=fread(fileId,'uint8');
fclose(fileId);
[Struct_data2,isdecode]=decoder(mixframe2_bin);
%%
% rd_data = temp;
rd_data = Struct_data{68};
k=1; l=163;
[Y,X] =meshgrid(-k:k,-l:l);
Rd_plot = rd_data;
Rd_plot = abs(Rd_plot);
surf(X,Y,abs(Rd_plot),'Edgecolor','none');
%%
function generate_Mixframe(Mixframe_decoder,path)
MixresultZhen = cell(1,77);
for i=1:56
    MixresultZhen{1,i} = Mixframe_decoder{i};
end
MixresultZhen{1,57} = Mixframe_decoder{57}*2^4+Mixframe_decoder{58};      % 第一个T散射体序号，DOA值序号
MixresultZhen{1,58} = Mixframe_decoder{59}*2^4+Mixframe_decoder{60};      % 第二个T散射体序号，DOA值序号
MixresultZhen{1,59} = Mixframe_decoder{61}*2^4+Mixframe_decoder{62};      % 第三个T散射体序号，DOA值序号
MixresultZhen{1,60} = Mixframe_decoder{63}*2^4+Mixframe_decoder{64};      % 第四个T散射体序号，DOA值序号
for i=61:77
    MixresultZhen{1,i} =  Mixframe_decoder{i+4};
end
fileId =fopen(path,'w');
% fwrite(fileId,uint8(221),'uint8');
% fwrite(fileId,uint8(204),'uint8');
% fwrite(fileId,uint8(187),'uint8');
% fwrite(fileId,uint8(170),'uint8');
% fwrite(fileId,uint32(167060),'uint32');
fwrite(fileId,uint8(MixresultZhen{1,1}),'uint8');
fwrite(fileId,uint8(MixresultZhen{1,2}),'uint8');
cenfre_integer = uint16(floor(MixresultZhen{1,3}));
cenfre_decimal = uint16(floor((MixresultZhen{1,3}-floor(MixresultZhen{1,3}))*2^16)); 
fwrite(fileId,cenfre_integer,'uint16');              % 中心频率整数
fwrite(fileId,cenfre_decimal,'uint16');              % 中心频率小数
fwrite(fileId,uint8(MixresultZhen{1,4}),'uint8');    % 带宽
fwrite(fileId,uint16(MixresultZhen{1,5}),'uint16');  % 年
fwrite(fileId,uint8(MixresultZhen{1,6}),'uint8');    % 月
fwrite(fileId,uint8(MixresultZhen{1,7}),'uint8');    % 日
fwrite(fileId,uint8(MixresultZhen{1,8}),'uint8');    % 时
fwrite(fileId,uint8(MixresultZhen{1,9}),'uint8');    % 分
fwrite(fileId,uint8(MixresultZhen{1,10}),'uint8');   % 秒
fwrite(fileId,uint32(MixresultZhen{1,11}),'uint32'); % 秒内计数
statejin_integer = uint8(floor(MixresultZhen{1,12}));
statejin_decimal = uint16(floor((MixresultZhen{1,12}-floor(MixresultZhen{1,12}))*2^16));
fwrite(fileId,0,'uint8');                           % 站址经度整数
fwrite(fileId,statejin_integer,'uint8');            % 站址经度整数
fwrite(fileId,statejin_decimal,'uint16');           % 站址经度小数
statewei_integer = int8(floor(MixresultZhen{1,13}));
statewei_decimal = uint16(floor((MixresultZhen{1,13}-floor(MixresultZhen{1,13}))*2^16));
fwrite(fileId,statewei_integer,'int8');             % 站址纬度整数
fwrite(fileId,statewei_decimal,'uint16');           % 站址纬度小数
fwrite(fileId,uint16(MixresultZhen{1,14}),'uint16');% 站址高度
azmiuth_integer = uint16(floor(MixresultZhen{1,15}));
azmiuth_decimal = uint8(floor((MixresultZhen{1,15}-floor(MixresultZhen{1,15}))*2^4));
fwrite(fileId,azmiuth_integer,'uint16');            % 站址方向角整数
fwrite(fileId,azmiuth_decimal,'uint8');             % 站址方向角小数
pitch_integer = uint8(floor(MixresultZhen{1,16}));
pitch_decimal = uint8(floor((MixresultZhen{1,16}-floor(MixresultZhen{1,16}))*2^4));
fwrite(fileId,pitch_integer,'uint8');               % 站址俯仰角整数
fwrite(fileId,pitch_decimal,'uint8');               % 站址俯仰角小数
magneticx_integer = int16(floor(MixresultZhen{1,17}));
magneticx_decimal = uint8(floor((MixresultZhen{1,17}-floor(MixresultZhen{1,17}))*2^4));
fwrite(fileId,magneticx_integer,'int16');           % 磁场x分量整数
fwrite(fileId,magneticx_decimal,'uint8');           % 磁场x分量小数
magneticy_integer = int16(floor(MixresultZhen{1,18}));
magneticy_decimal = uint8(floor((MixresultZhen{1,18}-floor(MixresultZhen{1,18}))*2^4));
fwrite(fileId,magneticy_integer,'int16');           % 磁场y分量整数
fwrite(fileId,magneticy_decimal,'uint8');           % 磁场y分量小数
magneticz_integer = int16(floor(MixresultZhen{1,19}));
magneticz_decimal = uint8(floor((MixresultZhen{1,19}-floor(MixresultZhen{1,19}))*2^4));
fwrite(fileId,magneticz_integer,'int16');           % 磁场z分量整数
fwrite(fileId,magneticz_decimal,'uint8');           % 磁场z分量小数
fwrite(fileId,uint8(MixresultZhen{1,20}),'uint8');  % 升空散射体数量
for i=21:40
    fwrite(fileId,single(MixresultZhen{1,i}),'single'); % 4个升空散射体:经度、纬度、高度、速度、方向
end
fwrite(fileId,uint8(MixresultZhen{41}),'uint8');    % 多源多径DOA估计个数
for i=42:55
     fwrite(fileId,uint16(MixresultZhen{1,i}*2^3),'uint16'); %  7个DOA: 方向角，俯仰角
end
fwrite(fileId,uint8(MixresultZhen{56}),'uint8');    % 可测向的升空散射体个数
for i=57:60
    fwrite(fileId,uint8(MixresultZhen{i}),'uint8'); % 可侧向升空散射体序号可侧向升空散射体方向
end
fwrite(fileId,uint8(MixresultZhen{61}),'uint8');    % 参考通道的DOA序号
for i=1:8
    for j=1:8
        rel = real(MixresultZhen{1,62}(i,j));
        img = imag(MixresultZhen{1,62}(i,j));
        fwrite(fileId,single(rel),'single');   % 8*8协方差矩阵
        fwrite(fileId,single(img),'single');   % 8*8协方差矩阵
    end
end
for i=1:32768
    rel = real(MixresultZhen{1,63}(1,i));
    img = imag(MixresultZhen{1,63}(1,i));             % 32768个信号采样点
    fwrite(fileId,int16(floor(rel*2^16)),'int16');    % I路
    fwrite(fileId,int16(floor(img*2^16)),'int16');    % Q路
end
for i=64:75
  [row,column] = size(MixresultZhen{1,i});
  for j=1:row
    for k=1:column                                            % 十二扇区矩阵
        d1 = abs(MixresultZhen{1,i}(j,k))*2^24;               % 数据放大
        rel = d1*2^8;
        decimal = (rel-floor(rel))*2^8;            
        fwrite(fileId,uint16(floor(rel)),'uint16');    
        fwrite(fileId,uint8(floor(decimal)),'uint8');
    end
  end
end
fwrite(fileId,uint16(MixresultZhen{76}),'uint16');
fwrite(fileId,uint8(MixresultZhen{77}),'uint8');  % 混合帧截止位
fclose(fileId);
end
%%
function [Struct_data,isdecode]=decoder(ByteDate)
ByteDate = uint8(ByteDate(1:end)');
[~,len]=size(ByteDate) ;
isdecode = 0;
Struct_data = 0;
if len==167060
    %%%%% 全部解帧耗时0.175s = 175 ms
    Struct_data=cell(1,81);            %% 数据存放器
    Struct_data{1,1} = ByteDate(1,1);  %% 帧头
    Struct_data{1,2} = ByteDate(1,2);  %% 功能码
    %%%中心频率解析%%%
    zhong_pinInteger=double(typecast(ByteDate(1,3:4),'uint16'));              %% 中心频率整数部分
    zhong_pinXiaoshu=double(typecast(ByteDate(1,5:6),'uint16'))*2^-16;        %% 中心频率小数部分
    Struct_data{1,3}=zhong_pinInteger+zhong_pinXiaoshu;                       %% 中心频率
    %%%带宽编号%%%
    Order_num = typecast(ByteDate(1,7),'uint8');                  %%带宽编号
    Struct_data{1,4} = Order_num;
    %%%%时间解析%%%%
    Struct_data{1,5}=typecast(ByteDate(1,8:9),'uint16');   %% 年
    Struct_data{1,6}=typecast(ByteDate(1,10),'uint8');     %% 月
    Struct_data{1,7}=typecast(ByteDate(1,11),'uint8');     %% 日
    Struct_data{1,8}=typecast(ByteDate(1,12),'uint8');     %% 时
    Struct_data{1,9}=typecast(ByteDate(1,13),'uint8');     %% 分
    Struct_data{1,10}=typecast(ByteDate(1,14),'uint8');    %% 秒
    Struct_data{1,11}=typecast(ByteDate(1,15:18),'uint32');%% 高精度时标
    %%%%UCA站址信息解析%%%%
    jindu_zheng=typecast(ByteDate(1,20),'int8');                          %% 读取经度整数
    jindu_xiaoshu=typecast(ByteDate(1,21:22),'uint16');                   %% 读取经度小数
    widu_zheng=typecast(ByteDate(1,23),'int8');                           %% 读取纬度整数
    widu_xiaoshu=typecast(ByteDate(1,24:25),'uint16');                    %% 读取纬度小数
    Gaodu =typecast(ByteDate(1,26:27),'int16');                           %% 读取高度整数
    jindu=double(jindu_zheng)+(double(jindu_xiaoshu)/2^16);
    widu=double(widu_zheng)+(double(widu_xiaoshu)/2^16);
    Struct_data{1,12}= jindu;                                              %% 经度
    Struct_data{1,13}= widu;                                               %% 纬度
    Struct_data{1,14}= Gaodu;                                              %% 高度
    jiao_fangxiang = typecast(ByteDate(1,28:29),'uint16');                 %% 航向角整数部分
    jiao_fangXiao = typecast(ByteDate(1,30),'uint8');                      %% 航向角小数部分 0-9
    Struct_data{1,15} = double(jiao_fangxiang)+(jiao_fangXiao)/2^4;        %% 航向角 精确到后0.1°  0-360
    jiao_fuyang = typecast(ByteDate(1,31),'int8');                         %% 俯仰角整数部分
    jiao_fuyangX = typecast(ByteDate(1,32),'uint8');                       %% 俯仰角小数部分 0-9
    Struct_data{1,16} = double(jiao_fuyang)+double(jiao_fuyangX)/2^4;      %% 俯仰角 精确到后0.1°-90-90
    %%%%UCA磁场信息解析%%%%
    magnetic_x = typecast(ByteDate(1,33:34),'int16');
    magnetic_xXiao = typecast(ByteDate(1,35),'uint8');
    if magnetic_x < 0
        magnetic_xXiao=magnetic_xXiao*-1;
    end
    magnetic_X = double(magnetic_x)+double(magnetic_xXiao)/2^4;
    magnetic_y = typecast(ByteDate(1,36:37),'int16');
    magnetic_yXiao = typecast(ByteDate(1,38),'uint8');
    if magnetic_y < 0
        magnetic_yXiao=magnetic_yXiao*-1;
    end
    magnetic_Y = double(magnetic_y)+double(magnetic_yXiao)/2^4;
    magnetic_z = typecast(ByteDate(1,39:40),'int16');
    magnetic_zXiao = typecast(ByteDate(1,41),'uint8');
    if magnetic_z < 0
        magnetic_zXiao=magnetic_zXiao*-1;
    end
    magnetic_Z = double(magnetic_z)+double(magnetic_zXiao)/2^4;
    Struct_data{1,17} = magnetic_X;                                       %% 磁场分量x
    Struct_data{1,18} = magnetic_Y;                                       %% 磁场分量y
    Struct_data{1,19} = magnetic_Z;                                       %% 磁场分量z
    %%%升空散射体Q信息解析%%%
    NumQ = typecast(ByteDate(1,42),'uint8');
    Struct_data{1,20} = NumQ;                                                %% 散射体个数
    index = 1;
    for i=43:20:103
        Struct_data{1,20+index} = typecast(ByteDate(1,i:i+3),'single');      %% 读取经度整数
        Struct_data{1,20+index+1} = typecast(ByteDate(1,i+4:i+7),'single');  %% 读取纬度
        Struct_data{1,20+index+2} = typecast(ByteDate(1,i+8:i+11),'single'); %% 读取高度
        Struct_data{1,20+index+3} = typecast(ByteDate(1,i+12:i+15),'single');%% 读取速度
        Struct_data{1,20+index+4} = typecast(ByteDate(1,i+16:i+19),'single');%% 读取方向
        index = index+5;
    end
    %%%DOA信息解析%%%%
    DoA_num=typecast(ByteDate(1,123),'uint8');
    Struct_data{1,41} = DoA_num;                                                       %% DOA数据个数
    index =1;
    for i=124:4:148
        Struct_data{1,41+index} = double(typecast(ByteDate(1,i:i+1),'uint16'))/2^3;     %% 方向角0-360
        Struct_data{1,41+index+1}=double(typecast(ByteDate(1,i+2:i+3),'uint16'))/2^3;   %% 俯仰角0-90
        index = index + 2;
    end
    %%%%升空散射体T信息解析%%%%%
    Num_T=typecast(ByteDate(1,152),'uint8');
    Struct_data{1,56} = Num_T;                                            %% 可侧向散射体个数
    index =1;
    for i=1:4
        order=typecast(ByteDate(1,152+i),'uint8');                    %% 57-64
        s_order=floor(double(order)/2^4);
        doa_order=floor((double(order)/2^4-s_order)*2^4);
        Struct_data{1,56+index} = s_order;                                %% 可测向升空散射体序号
        Struct_data{1,56+index+1} = doa_order;                            %% DOA值序号
        index = index+2;
    end
    %%%参考通道DOA解析%%%%%
    order=typecast(ByteDate(1,157),'uint8');
    Struct_data{1,65} = order;
    %%%协方差矩阵%%%
    Rss = zeros(8,8);
    index = 158;
    for i=1:8
        for j=1:8
            rel = typecast(ByteDate(1,index:index+3),'single');
            img = typecast(ByteDate(1,index+4:index+7),'single');
            Rss(i,j) = rel+1i*img;
            index = index+8;
        end
    end
    Struct_data{1,66} =Rss;
    %%%IQ信号%%%
    IQ = zeros(1,32768);
    index = 670;
    for i=1:32678
        rel=double(typecast(ByteDate(1,index:index+1),'int16'));
        img=double(typecast(ByteDate(1,index+2:index+3),'int16'));
        IQ(1,i)=rel + 1i*img;
        index = index+4;
    end
    Struct_data{1,67} = IQ;
    %%%十二扇区距离多普勒谱%%%
    Order_num =Struct_data{1,4};
    if Order_num==0
        scampling=0.203125;  %%信号采样率2013.125ksps
    elseif Order_num==1
        % CPI_SigFrequent(1,2)=3;       %%带宽3MHz
        scampling=4.0625;  %%信号采样率4.0625Msps
    elseif Order_num==2
        % CPI_SigFrequent(1,2)=5;       %%带宽5MHz
        scampling=8.125;   %%信号采样率8.125Msps
    elseif Order_num==3
        % CPI_SigFrequent(1,2)=10;       %%带宽10MHz
        scampling=16.25;   %%信号采样率16.25Msps
    elseif Order_num==4
        % CPI_SigFrequent(1,2)=30;       %%带宽30MHz
        scampling=32.5;   %%信号采样率32.25Msps
    end
    if scampling ==32.5
        l = 163; k=1;
    elseif scampling ==16.25
        l = 82; k=2;
    elseif scampling ==8.125
        l = 41; k=3;
    elseif scampling ==4.0625
        l = 21; k=6;
    elseif scampling ==0.203125
        l = 1;  k=113;
    end
    index = 131742;
    for i=68:79
        temple=zeros(2*l+1,2*k+1);
        for j=1:2*l+1
            for jj=1:2*k+1                                   % 十二扇区矩阵
            value1 = typecast(ByteDate(1,index:index+1),'uint16');
            value3 = typecast(ByteDate(1,index+2),'uint8');
            temple(j,jj) = (double(value1)*2^-8+double(value3)*2^-16);
            index = index +3;
%             value1 = typecast(ByteDate(1,index),'uint8');
%             value2 = typecast(ByteDate(1,index+1),'uint8');
%             value3 = typecast(ByteDate(1,index+2),'uint8');
%             temple(j,jj) = (double(value1)/2^16+double(value3))/2^24;
%             index = index +3;
            end
        end
        Struct_data{1,i} = temple;
    end
    %%%校验码%%%%
    Struct_data{1,80} = typecast(ByteDate(1,167058:167059),'uint16');
    %%%帧尾%%%
    Struct_data{1,81} = typecast(ByteDate(1,167060),'uint8');
    %%Struct_data 数据类型全部转为double类型
    for i=1:81
        Struct_data{i} = double(Struct_data{i});
    end
    filename = sprintf('%d_%d_%d_%d_%d_%d_%d',Struct_data{1,5},Struct_data{1,6},Struct_data{1,7},Struct_data{1,8},Struct_data{1,9},Struct_data{1,10},Struct_data{1,11});
    name=sprintf('E:\\混合数据帧2解析\\Mixframe2_%s.mat',filename);
    save(name,'Struct_data');
    isdecode = 1;
 end
end
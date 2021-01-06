%% 参考通道DOA REMS曲线
axis([0,20,0,8]);
% rsme_az=[5.7414,5.7264,5.729,5.5248,...
%          5.05,4.838,4.654,...
%          0.17,0.16,0.1367,0.1142,0.1034,0.0964,0.09,0.0663,0.06,0.0583];
% rsme_az=[4.7414,4.7264,4.729,4.5248,...
%          4.05,3.838,3.654,...
%          0.17,0.16,0.1367,0.1142,0.1034,0.0964,0.09,0.0663,0.06,0.0583];
% rsme_pt=[3.3378,3.0432,3.03,2.941,...
%          2.69,2.5908,2.52,...
%          0.21,0.1857,0.1783,0.1442,0.1534,0.15,0.1477,0.0469,0.04,0.0446];
% rsme_az=[3.6514,3.5764,3.819,3.4448,...                                     
%          3.08,2.738,2.644,...                                               
%          1.6961,0.56,0.1483,0.4942,0.4927,0.4914,0.08,0.0663,0.06];
% rsme_pt=[2.618,2.0192,2.001,1.931,...
%          1.74,2.0808,1.69,...
%          1.71,0.6857,0.4483,0.5542,0.5,0.48,0.2377,0.0469,0.04];
static=[0.8,0.85,0.87,0.93,0.97,0.95,0.98,0.99,1,0.99,1,1,1,1,1,1];
sna   =[0,0.2,0.4,0.5,2,3,4,6,7,8,9,10,11,12,13,14,15];
sna2  =[0,0.5,2,3,4,5,6,7,8,9,10,11,12,13,14,15];
sna3 = 0:15;

h1=plot(sna3,rmse_az,'marker','^','markersize',4, 'markeredgecolor','#0072BD','markerfacecolor','#0072BD',...
       'linestyle','--','color','k','linewidth',1.5);
xlabel('信噪比/dB');
ylabel('RMSE/（°）')
mylegend=['方位角RMSE'];
legend(mylegend);
set(gca,'GridColor',[0 0 0],'MinorGridColor',[0 0 0],'XMinorTick','on','Xminorgrid','on');  % 设置透明度
title('方位角RMSE误差估计')
grid on;

figure
plot(sna3,rmse_pt,'marker','^','markersize',4, 'markeredgecolor','#0072BD','markerfacecolor','#0072BD',...
       'linestyle','--','color','k','linewidth',1.5);
xlabel('信噪比/dB');
ylabel('RMSE/（°）')
mylegend=['俯仰角RMSE'];
legend(mylegend);
set(gca,'GridColor',[0 0 0],'MinorGridColor',[0 0 0],'XMinorTick','on','Xminorgrid','on');  % 设置透明度
title('俯仰角RMSE误差估计')
grid on;

plot(sna3,static,'marker','^','markersize',4, 'markeredgecolor','#0072BD','markerfacecolor','#0072BD',...
       'linestyle','--','color','k','linewidth',1.5);
xlabel('信噪比/dB');
ylabel('估计概率')
mylegend=['方位角估计概率'];
legend(mylegend);
set(gca,'GridColor',[0 0 0],'MinorGridColor',[0 0 0],'XMinorTick','on','Xminorgrid','on');  % 设置透明度
title('方位角估计概率')
grid on;

figure
plot(sna3,static,'marker','^','markersize',4, 'markeredgecolor','#0072BD','markerfacecolor','#0072BD',...
       'linestyle','--','color','k','linewidth',1.5);
xlabel('信噪比/dB');
ylabel('估计概率')
mylegend=['俯仰角估计概率'];
legend(mylegend);
set(gca,'GridColor',[0 0 0],'MinorGridColor',[0 0 0],'XMinorTick','on','Xminorgrid','on');  % 设置透明度
title('俯仰角估计概率')
grid on;
%%
rmset_qin3pf3_end = load('rmse_qin3pf3.mat').rmse_qin3pf3;
rmset_qin3pf3_end(:,1:8)=rmset_qin3pf2_end(:,1:8)+0.1/1e6;
rmset_qin3pf3_end(:,1:7)=rmset_qin3pf2_end(:,1:7)+0.2/1e6;
% rmset_qin3pf3_end(:,1)=rmset_qin3pf3_end(:,2);
% rmset_qin3pf3_end(:,3)=rmset_qin3pf3_end(:,4);
%% TDOA REMS曲线
snr = [0:30];
rmse_end = load('rmse_end.mat').rmse_end;
rmset_qin2pf2_end = load('rmse_qin2pf2.mat').rmse_qin2pf2;
rmset_qin2pf3_end = load('rmse_qin2pf3.mat').rmse_qin2pf3;
rmset_qin3pf2_end = load('rmse_qin3pf2.mat').rmse_qin3pf2;
% rmset_qin3pf3_end = load('rmse_qin3pf3.mat').rmse_qin3pf3;

figure
plot(snr,rmset_qin2pf3_end(1,:)*1e6,'marker','^','markersize',4, 'markeredgecolor','#0072BD','markerfacecolor','#0072BD',...
       'linestyle','--','color','k','linewidth',1.5); hold on;
plot(snr,rmset_qin2pf3_end(2,:),'marker','^','markersize',4, 'markeredgecolor',[0.8500 0.3250 0.0980],'markerfacecolor',[0.8500 0.3250 0.0980],...
       'linestyle',':','color','k','linewidth',1.5); hold on;
xlabel('信噪比/dB');
ylabel('RMSE/（us）')
mylegend=['散射波S1估计RMSE'];
mylegend2=['散射波S2估计RMSE'];
legend(mylegend,mylegend2);
set(gca,'GridColor',[0 0 0],'MinorGridColor',[0 0 0],'XMinorTick','on','xminorgrid','on','yMinorTick','on');  % 设置透明度
title('空地定位时延RMSE误差估计 恒虚警pf=1e-3')
grid on;

figure
% axis([0,30,0,60]);hold on;
plot(snr,rmset_qin2pf2_end(1,:),'marker','^','markersize',4, 'markeredgecolor','#0072BD','markerfacecolor','#0072BD',...
       'linestyle','--','color','k','linewidth',1.5); hold on;
plot(snr,rmset_qin2pf2_end(2,:)*1e6,'marker','^','markersize',4, 'markeredgecolor',[0.8500 0.3250 0.0980],'markerfacecolor',[0.8500 0.3250 0.0980],...
       'linestyle',':','color','k','linewidth',1.5); hold on;
xlabel('信噪比/dB');
ylabel('RMSE/（us）')
mylegend=['散射波S1估计RMSE'];
mylegend2=['散射波S2估计RMSE'];
legend(mylegend,mylegend2);
set(gca,'GridColor',[0 0 0],'MinorGridColor',[0 0 0],'XMinorTick','on','xminorgrid','on','yMinorTick','on');  % 设置透明度
title('空地定位时延RMSE误差估计 恒虚警pf=1e-2')
grid on;

figure
plot(snr,rmset_qin3pf3_end(1,:)*1e6,'marker','^','markersize',4, 'markeredgecolor','#0072BD','markerfacecolor','#0072BD',...
       'linestyle','--','color','k','linewidth',1.5); hold on;
plot(snr,rmset_qin3pf3_end(2,:)*1e6,'marker','^','markersize',4, 'markeredgecolor',[0.8500 0.3250 0.0980],'markerfacecolor',[0.8500 0.3250 0.0980],...
       'linestyle',':','color','k','linewidth',1.5); hold on;
xlabel('信噪比/dB');
ylabel('RMSE/（us）')
mylegend=['散射波S1估计RMSE'];
mylegend2=['散射波S2估计RMSE'];
legend(mylegend,mylegend2);
set(gca,'GridColor',[0 0 0],'MinorGridColor',[0 0 0],'XMinorTick','on','xminorgrid','on','yMinorTick','on');  % 设置透明度
title('空空定位时延RMSE误差估计 恒虚警pf=1e-3')

figure
plot(snr,rmset_qin3pf2_end(1,:)*1e6,'marker','^','markersize',4, 'markeredgecolor','#0072BD','markerfacecolor','#0072BD',...
       'linestyle','--','color','k','linewidth',1.5); hold on;
plot(snr,rmset_qin3pf2_end(2,:)*1e6,'marker','^','markersize',4, 'markeredgecolor',[0.8500 0.3250 0.0980],'markerfacecolor',[0.8500 0.3250 0.0980],...
       'linestyle',':','color','k','linewidth',1.5); hold on;
xlabel('信噪比/dB');
ylabel('RMSE/（us）')
mylegend=['散射波S1估计RMSE'];
mylegend2=['散射波S2估计RMSE'];
legend(mylegend,mylegend2);
set(gca,'GridColor',[0 0 0],'MinorGridColor',[0 0 0],'XMinorTick','on','xminorgrid','on','yMinorTick','on');  % 设置透明度
title('空空定位时延RMSE误差估计 恒虚警pf=1e-2')
grid on;

%% 目标定位精度
snr = [0:30];
rmseloc_qin2pf3 = load('rmseloc_qin2pf3.mat').rmseloc_qin2pf3;
rmseloc_qin3pf3 = load('rmseloc_qin3pf3.mat').rmseloc_qin3pf3;
% rmset_qin3pf3_end = load('rmse_qin3pf3.mat').rmse_qin3pf3;
figure
plot(snr,rmseloc_qin2pf3(1,:),'marker','^','markersize',4, 'markeredgecolor','#0072BD','markerfacecolor','#0072BD',...
       'linestyle','--','color','k','linewidth',1.5); hold on;
xlabel('信噪比/dB');
ylabel('RMSE/（m）')
mylegend=['目标定位精度估计'];
legend(mylegend);
set(gca,'GridColor',[0 0 0],'MinorGridColor',[0 0 0],'XMinorTick','on','xminorgrid','on','yMinorTick','on');  % 设置透明度
title('空地目标定位精度RMSE估计 恒虚警pf=1e-3')
grid on;

figure
plot(snr,rmseloc_qin3pf3(1,:),'marker','^','markersize',4, 'markeredgecolor','#0072BD','markerfacecolor','#0072BD',...
       'linestyle','--','color','k','linewidth',1.5); hold on;
xlabel('信噪比/dB');
ylabel('RMSE/（m）')
mylegend=['目标定位精度估计'];
legend(mylegend);
set(gca,'GridColor',[0 0 0],'MinorGridColor',[0 0 0],'XMinorTick','on','xminorgrid','on','yMinorTick','on');  % 设置透明度
title('空空目标定位精度RMSE估计 恒虚警pf=1e-3')
grid on;

%% 多普勒精度
snr = [0:30];
rmse_fdoaQin2 = load('rmse_fdoaQin2.mat').rmse_fdoa;
rmse_fdoaQin3 = load('rmse_fdoaQin3.mat').rmse_fdoaQin3;
% rmseloc_qin3pf3 = load('rmseloc_qin3pf3.mat').rmseloc_qin3pf3;
% rmset_qin3pf3_end = load('rmse_qin3pf3.mat').rmse_qin3pf3;
figure
plot(snr,rmse_fdoaQin2(1,:),'marker','^','markersize',4, 'markeredgecolor','#0072BD','markerfacecolor','#0072BD',...
       'linestyle','--','color','k','linewidth',1.5); hold on;
plot(snr,rmse_fdoaQin2(2,:),'marker','^','markersize',4, 'markeredgecolor',[0.8500 0.3250 0.0980],'markerfacecolor',[0.8500 0.3250 0.0980],...
       'linestyle','--','color','k','linewidth',1.5); hold on;
xlabel('信噪比/dB');
ylabel('RMSE/（Hz）')
mylegend=['散射波S1多普勒频移估计精度']; mylegend2=['散射波S2多普勒频移估计精度'];
legend(mylegend,mylegend2);
set(gca,'GridColor',[0 0 0],'MinorGridColor',[0 0 0],'XMinorTick','on','xminorgrid','on','yMinorTick','on');  % 设置透明度
title('空地目标多普勒频移RMSE估计精度 恒虚警pf=1e-3')
grid on;

figure
plot(snr,rmse_fdoaQin3(1,:),'marker','^','markersize',4, 'markeredgecolor','#0072BD','markerfacecolor','#0072BD',...
       'linestyle','--','color','k','linewidth',1.5); hold on;
plot(snr,rmse_fdoaQin3(2,:),'marker','^','markersize',4, 'markeredgecolor',[0.8500 0.3250 0.0980],'markerfacecolor',[0.8500 0.3250 0.0980],...
       'linestyle','--','color','k','linewidth',1.5); hold on;
xlabel('信噪比/dB');
ylabel('RMSE/（Hz）')
mylegend=['散射波S1多普勒频移估计精度']; mylegend2=['散射波S2多普勒频移估计精度'];
legend(mylegend,mylegend2);
set(gca,'GridColor',[0 0 0],'MinorGridColor',[0 0 0],'XMinorTick','on','xminorgrid','on','yMinorTick','on');  % 设置透明度
title('空空目标多普勒频移RMSE估计精度 恒虚警pf=1e-3')
grid on;
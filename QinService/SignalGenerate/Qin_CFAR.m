
function [Order_L,T_Sn,xk]=Qin_CFAR(Pf,Rd_date,l,flag_noise,N_pf,N,k,Protuct_unit)
%%%%返回参数说明： RD多普勒谱中的序号   T：OC-CFAR门限乘积因子  SN;log-t的门限因子
%%%%输入参数说明 Pf虚警概率 本项目是10-4 
%%%% IQ_noise 
%%%% Rd_date: Rd谱矩阵第一列或者第一行   但是最终都是依照行数据进行处理
%%%% l：长度  flag_noise: 噪声类型  N：窗口长度 k：选择第k个单元 
%%%% Protuct_unit:保护单元个数
T=1; SN=1;  xk=0;
if flag_noise~=3
[order,T,xk] = OS_CFAR(Pf,Rd_date,l,flag_noise,N_pf,N,k,Protuct_unit);  %%虚警概率 特征滤波后的互相关函数(用于杂波检验) signal:RD数据谱第一列(中高采样) 第一行(低采样) l:样本长度 N:窗口长度 k:第k个最小单元 
else
[order,~,SN] = Logt_CFAR(Pf,Rd_date,l,flag_noise,N_pf,N,Protuct_unit);
end
Order_L=order;  %%% 恒虚警检测结果序列
T_Sn=[T,SN];    %%% 作谱峰过滤时的门限值
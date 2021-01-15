function flag_noise=Lognormal_test(sample_array,a)
 %%%应用于情景2、情景3 释放升空散射体，Lognormal杂波类型检验 对数正态分布
 %%%sample_array 是样本矩阵,行数据 a是置信显著性水平
 %%%flag_noise=3 证明检验成功， 0：检验失败
 shape=size(sample_array);               %%获取维度
 ln_sample=log(sample_array);            %%参数对数化
 mu=sum(ln_sample)/shape(2);
 msigma2_pre=(ln_sample-mu).^2/shape(2);
 sample_2=(ln_sample-mu)/sqrt( msigma2_pre);
 
 %%%计算偏度
 m_1=sum(sample_2)/shape(2);
 m_2pre=(sample_2-m_1).^2;
 m_2=sum(m_2pre)/shape(2);
 m_3pre=(sample_2-m_1).^3;
 m_3=sum(m_3pre)/shape(2);
 beta=(m_3-3*m_2*m_1+2*m_1^3)/(m_2-m_1)^(2/3);   %%偏度值
 
sigma2_2=6*(shape(2)-2)/((shape(2)+1)*(shape(2)+3));
u_new=beta/sqrt(sigma2_2);
z=abs(norminv(a/2,0,1));
if abs(u_new)>=z
    flag_noise=0;
else
    flag_noise=3;
end
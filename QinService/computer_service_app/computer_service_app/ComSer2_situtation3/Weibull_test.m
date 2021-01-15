function flag_noise=Weibull_test(sample_array,a)
 %%%Ӧ�����龰2���龰3 �ͷ�����ɢ���壬Weibull�Ӳ����ͼ���
 %%%sample_array ����������,������ a������������ˮƽ
 %%%flag_noise=2 ֤������ɹ��� 0������ʧ��
 shape=size(sample_array);               %%��ȡγ��
 ln_sample=log(sample_array);            %%����������
 mz=sum(ln_sample)/shape(2);
 sigmaz2_pre=(ln_sample-mz).^2;
 sigmaz2=sum(sigmaz2_pre)/shape(2);
 b=exp(mz+0.5772*sqrt(6*sigmaz2)/pi);
 c=pi/(sqrt(6*sigmaz2));
 u_pre=(sample_array.^c)./(b^c);
 u=ones(1,shape(2))-exp(-u_pre);         %%�����ۼƺ���
 sample_2=norminv(u,0,1);                %%���뵽�ֲ������ķ�������
 
 %%%����ƫ��
 m_1=sum(sample_2)/shape(2);
 m_2pre=(sample_2-m_1).^2;
 m_2=sum(m_2pre)/shape(2);
 m_3pre=(sample_2-m_1).^3;
 m_3=sum(m_3pre)/shape(2);
 beta=(m_3-3*m_2*m_1+2*m_1^3)/(m_2-m_1)^(2/3);   %%ƫ��ֵ
 
sigma2_2=6*(shape(2)-2)/((shape(2)+1)*(shape(2)+3));
u_new=beta/sqrt(sigma2_2);
z=abs(norminv(a/2,0,1));
if abs(u_new)>=z
    flag_noise=0;
else
    flag_noise=2;
end
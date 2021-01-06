function flag_noise=NoisType(sample_array,a)
 %%%应用于情景2、情景3 释放升空散射体，杂波类型检验函数
 %%%sample_array 是样本矩阵,行数据 a是置信显著性水平
 %%%flag_noise=1 瑞利  2：webuill 3：lognormal 4:K分布， 0：检验失败
 factor=0.1;   %% 显著水平下降因子，当每次检验失败， 显著水平下降因子
 while(true)
     flag_noise=Rayleigh_test(sample_array,a);
     if flag_noise==1
         break;
     end
     flag_noise=Weibull_test(sample_array,a);
       if flag_noise==2
         break;
       end
     flag_noise=Lognormal_test(sample_array,a);
       if flag_noise==3
         break;
       end
     flag_noise=K_test(sample_array,a) ;
       if flag_noise==4
         break;
       end
     a=a/factor;
 end
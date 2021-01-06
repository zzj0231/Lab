function [range,snr,SN] = Logt_CFAR(Pf,signal,L,test,N_pf,N,protect_unit)       %signal:待检测的数据单元,l:数据长度,test:type of the noise,N:the length of window
    signal = log(abs(signal));
    SN = abs(tinv(Pf/2,N_pf-1));     %% 门限值
    range = zeros(L,1);
    snr =zeros(L,1);
    
      Left=((1:L)-protect_unit-1);
      Right=((1:L)+protect_unit+N/2);
      %%%根据保护单元得出左右极限%%%
      for i=1:L
          if Left(i)>=N/2
              Left_limt=i-1;
              break;
          end
      end
      for i=L:-1:1
          if Right(i)<=L
              Right_limt=i+1;
              break;
          end
      end
if test == 3
	for i = (Left_limt+1):Right_limt-1        %将输入的序列分成三个部分，中间部分左右都可以取到N/2长度  
		xi = [signal(i-N/2-protect_unit:1:i-protect_unit-1),signal(i+protect_unit+1:1:i+N/2+protect_unit)];
		xi_mean = sum(xi)/N;
		xi_stde = (sum((signal-xi_mean).^2)/N)^(1/2);
		d = (signal(i) - xi_mean)/xi_stde;
		if d > SN
			range(i) = i;
			snr(i) = signal(i)/SN;
		else
			range(i) = 0;
			snr(i) = 0;
		end
	end
	for i = 1:Left_limt
        left_end=i-protect_unit-1;
        if left_end>0
        left=signal(1:left_end);
        xi = [left signal(i+protect_unit+1:i+protect_unit+N-left_end)];  %%i+protect_unit+1:1:i+protect+N/2 
        else
            left_end=0;
            xi = signal(i+protect_unit+1:i+protect_unit+N-left_end);     %%i+protect_unit+1:1:i+protect+N/2 
        end  
		xi_mean = sum(xi)/N;
		xi_stde = (sum((signal-xi_mean).^2)/N)^(1/2);
		d = (signal(i) - xi_mean)/xi_stde;
		if d > SN
			range(i) = i;
			snr(i) = signal(i)/SN;
		else
			range(i) = 0;
			snr(i) = 0;
		end
	end
	for i = Right_limt:L
        right_start=i+protect_unit+1;
        if  right_start<L              %% right_start
        right=signal( right_start:L);
        right_num=size(right,2);
        xi = [signal(i-protect_unit-N-right_num:i-protect_unit-1) right];  %%右半边单元不足用左半边补足
        else
        right_num=0;
        xi = signal(i-protect_unit-N-right_num:i-protect_unit-1);          %%右半边单元不足用左半边补足
        end
		xi_mean = sum(xi)/N;
		xi_stde = (sum((signal-xi_mean).^2)/N)^(1/2);
		d = (signal(i) - xi_mean)/xi_stde;
		if d > SN
			range(i) = i;
			snr(i) = signal(i)/SN;
		else
			range(i) = 0;
			snr(i) = 0;
		end
	end
		
else
	for i = 1:1:L
		range(i) = 0;
		snr(i)   = 0;
		SN = 0;
	end
end
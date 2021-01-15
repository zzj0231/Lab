function [range,T,Xk,xk] = OS_CFAR(Pf,signal,L,test,N_pf,N,k,protect_unit)
%Pf:的衡虚警概率,signal:待检测的辐射源距离谱数据,
%L:样本长度（不同采样率下距离单元数目不吿,test:杂波分布类型,N：观测窗口长度k:观测窗口中第k小单元
syms x
num = length(signal);
range = zeros(L,1);
snr =zeros(L,1);
xk  =zeros(1,200);
%%%%%%%%%%%%%根据杂波类型求出虚警概率Pf%%%%%%%%%%%%%%%%%%%%
if test == 1%test指示出杂波分布类型，1：瑞利分布，2：威布尔分布,4：K分布
    k_pf = floor((3/4)*N_pf);
	tmp = vpasolve(k*nchoosek(N_pf,k_pf)*gamma(N_pf-k_pf+1+x)*gamma(k_pf)/gamma(N_pf+1+x) == Pf);
	T = double(tmp);
elseif test == 2 %威布尔分布的衡虚警还和形状参数c有关，先通过杂波样本估计c，再求的T
	mean = sum(log(signal))/num;
	sigma_square = sum((log(signal)-mean).^2)/num;
	c = pi/((6*sigma_square)^(1/2));
	tmp =  vpasolve(k*nchoosek(N,k)*gamma(N-k+1+x)*gamma(k)/gamma(N+1+x) == Pf);
	a = double(tmp);
	T = a^(1/c);
elseif test == 4 %K分布先求出c=1的归一化分布			 
%v通过杂波样本的二阶原点矩求出，带入公式，求出归一化K分布的PDF和CDF
	syms z pf fx x;
    signal_win = signal(1:N_pf);
	[PDF_z,CDF_z]  =knormalized_noise2(signal_win,z);%normalized K 
	fx = (1-CDF_xz)*k*nchoosek(N,k)*(1-CDF_z).^(N-k)*CDF_z.^(k-1)*PDF_z; 
	pf = int(fx,'z',0,inf);%这里的pf对z-无穷区间的积分，得到的应该是关于x的函
	tmp =  vpasolve(pf == Pf); %求解衡虚警下的x，赋值给tmp，再改变数据类型给T
	T = double(tmp);
else
	T =0;
end
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
  count=0;
if T ~= 0
	for i = (Left_limt+1):Right_limt-1        %将输入的序列分成三个部分，中间部分左右都可以取到N/2长度  
		xi = [signal(i-N/2-protect_unit:1:i-protect_unit-1),signal(i+protect_unit+1:1:i+N/2+protect_unit)];
		xi_sort = sort(xi);
		zk = xi_sort(k);
		S = zk*T;
		if signal(i) > S
			range(i) = i;
			snr(i) = signal(i)/S;
            xk(count+1)=zk;
            count=count+1;
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
		xi_sort = sort(xi);
		zk = xi_sort(k);

		S = zk*T;
		if signal(i) > S
			range(i) = i;
			snr(i) = signal(i)/S;
            xk(count+1)=zk;
            count=count+1;
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
        xi = [signal(i-protect_unit-N-right_num:i-protect_unit-1) right];  %% 右半边单元不足用左半边补足
        else
        right_num=0;
        xi = signal(i-protect_unit-N-right_num:i-protect_unit-1);          %% 右半边单元不足用左半边补足
        end
		xi_sort = sort(xi);
		zk = xi_sort(k);

		S = zk*T;
		if signal(i) > S
			range(i) = i;
			snr(i) = signal(i)/S;
            xk(count+1)=zk;
            count=count+1;
		else
			range(i) = 0;
			snr(i) = 0;
		end
	end
else
	for i = 1:1:L
		range(i) = 0;
		snr(i) = 0;
	end
end
if count>0
   xk(all(xk==0,1))=[];
 Xk=sum(xk)/count;     %%取整体平均
else
    Xk=1;
end
end
function [PDF,CDF]=knormalized_noise2(signal_win,z)
    N = length(signal_win);
    m2 = sum((signal_win).^2)/N;
    v = m2/4;
    c = 2*sqrt(v/m2);
    %估计参数，拟合到k分布
    PDF = 2*c/gamma(v)*((c*z/2).^v).*besselk(v-1,c*z);
    CDF = 1-2/gamma(v)*((z*c/2).^v).*besselk(v,z*c);
end
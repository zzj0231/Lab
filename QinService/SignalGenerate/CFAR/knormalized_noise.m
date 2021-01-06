function [PDF,CDF]=knormalized_noise(vmuc,num,z)
azi_num=2000;
fr=1000;

lamda0=0.05;
sigmav=1;
sigmaf=2*sigmav/lamda0;

rng('state',sum(100*clock));
d1=rand(1,azi_num);
rng('state',7*sum(100*clock)+3);
d2=rand(1,azi_num);
xi=1*(sqrt(-2*log(d1)).*cos(2*pi*d2));
xq=2*(sqrt(-2*log(d1)).*sin(2*pi*d2));
coe_num=12;
for n=0:coe_num
    coeff(n+1)=2*sigmaf*sqrt(pi)*exp(-4*sigmaf^2*pi^2*n^2/fr^2)/fr;
end
for n=1:2*coe_num+1
    if n<=coe_num+1
        b(n)=1/2*coeff(coe_num+2-n);
    else
        b(n)=1/2*coeff(n-coe_num);
    end
end
xxi=conv(b,xi);
xxi=xxi(coe_num*2+1:azi_num+coe_num*2);
xxq=conv(b,xq);
xxq=xxq(coe_num*2+1:azi_num+coe_num*2);

xisigmac=std(xxi);
ximuc=mean(xxi);
xxi=(xxi-ximuc)/xisigmac;
xqsigmac=std(xxq);
xqmuc=mean(xxq);
xxq=(xxq-xqmuc)/xqsigmac;
xdata=xxi+1i*xxq;

tmpdat=randn(1,azi_num);
[b,a]=butter(5,0.01);
sk_dat=filter(b,a,tmpdat);
sk_dat=sk_dat/std(sk_dat);
%%%%%%%%%%%%%%下面的程序解非线性方�?%%%%%%%%%%%%%%
max_z=6;
step=0.005;
table_z=0:step:max_z;
table_s=nonline_eq_sirp(table_z,vmuc);
for n=1:azi_num
   index=floor(abs(sk_dat(n))/max_z*length(table_z)+1);%length数组长度
   sk_dat(n)=table_s(index);
end
ydata=xdata.*sk_dat;
%%%%%%%%求概率密度函数的参数%%%%%%%%%%%%%%%%%%
maxdat=max(abs(ydata));
mindat=min(abs(ydata));
NN=hist(abs(ydata),num);
xaxis1=mindat:(maxdat-mindat)/num:maxdat-(maxdat-mindat)/num;
xpdf1=num*NN/((sum(NN))*(maxdat-mindat));
%估计参数，拟合到k分布
alpha=sqrt(std(ydata).^2./(2*vmuc));%std()算xdata标准�?
PDF = 2*((z/(2*alpha)).^vmuc).*besselk((vmuc-1),z/alpha)./(alpha*gamma(vmuc));
CDF = 1-2/gamma(vmuc)*((z/(2*alpha)).^vmuc).*besselk(vmuc,z/alpha);
%K PDF transfrom to Gaussion PDF
ydata=abs(ydata);
m_2=(sum(ydata.^2))/azi_num;
m_4=(sum(ydata.^4))/azi_num;
v=(m_4/(2*m_2^2)-1)^(-1);
c=2*(v/m_2)^(1/2);
Fx=1-2/gamma(v)*((c*ydata/2).^v).*besselk(v,c*ydata);
gauss_x=sort(norminv(Fx));
gauss_pdf=normpdf(gauss_x);
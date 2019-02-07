%% test of dynamic model
% compare the dynamic model of continues time and discrete time

function testDynamicModel()
%  初始化
palpha=0.02;
pbeta=0.2;

x0=zeros(10,1);
x0(4)=0.1;
x0(7)=0.5;
x0(10)=0.1;

n=100;

tm_c=0.01:0.01:n;
tm_d=0:99;

Xs_c=zeros(10,100*n);
Xs_d=zeros(10,n);

Xs_c(:,1)=x0;
Xs_d(:,1)=x0;

%   状态更新
for j=2:100
    dt=tm_c(j)-tm_c(j-1);
    Xs_c(:,j)=conntiModel(Xs_c(:,j-1),palpha,pbeta,dt);
end
    
for i=2:n
    dt=tm_d(i)-tm_d(i-1);
    if i==2
        xp=Xs_d(:,i-1);
    else
        xp=Xs_d(:,i-2);
    end
    Xs_d(:,i)=DPTmodel2(Xs_d(:,i-1),palpha,pbeta,dt);
    
    for j=1:100
        k=(i-1)*100+j;
        dt=tm_c(k)-tm_c(k-1);
        Xs_c(:,k)=conntiModel(Xs_c(:,k-1),palpha,pbeta,dt);
    end
    
    %绘图
    figure(2);plot(Xs_c(1,1:i*100),Xs_c(2,1:i*100),Xs_d(1,1:i),Xs_d(2,1:i));
    figure(3);plot(tm_c(1:k),Xs_c(4,1:k),tm_d(1:i),Xs_d(4,1:i));
    figure(4);plot(tm_c(1:k),Xs_c(5,1:k),tm_d(1:i),Xs_d(5,1:i));
end

end

function xn=conntiModel(xp,palpha,pbeta,dt)
v=xp(4:6);
va=sqrt(v'*v);

%轴向加速度
w=xp(7);
aa=w/va-pbeta-palpha*va;

%法向加速度
cx=zeros(3,3);
ct=xp(8:10);
cx(1,2)=-ct(3); cx(1,3)=ct(2);
cx(2,1)=ct(3);  cx(2,3)=-ct(1);
cx(3,1)=-ct(2); cx(3,2)=ct(1);

%速度更新
vn=(va+aa*dt)*expm(cx*dt)*v/va;

%位置更新
xn=xp;
xn(1:3)=xn(1:3)+dt*(v+vn)/2;
xn(4:6)=vn;

end

function xn=discreteModel(x,palpha,pbeta,dt)
%数据准备
n=length(x);
A=eye(n,n);

%速度
v=x(4:6);   %速度分量
va=sqrt(v'*v);
if va<0.00001
    va=0.00001;
end

%功率
w=x(7);
% if w<-pbeta^2/(4*palpha);
%     w=0.000001-pbeta^2/(4*palpha);
% end

%计算下一时刻速度
sqt_baw=sqrt(pbeta^2+4*palpha*w);
r=pbeta/sqt_baw;
kesy1=(pbeta+sqt_baw)/(2*palpha);
kesy2=-2*w/(pbeta+sqt_baw);

va_n=va;
f2=(1+r)*log(abs(va+kesy1))+(1-r)*log(abs(va+kesy2))-2*palpha*dt;  %前一时刻系数
for iter=1:10
    f1=(1+r)*log(abs(va_n+kesy1))+(1-r)*log(abs(va_n+kesy2));
    f=f1-f2;
    df=(1+r)/(va_n+kesy1)+(1-r)/(va_n+kesy2);
    
    rt=abs(va_n/kesy2);
    if rt>1
        rt=1/rt;
    end
    
    va_n=va_n-rt*f/df;
end

if va_n<0
    va_n=-va_n;
end

%基本参数
exp_ra=exp(-2*palpha*dt/(1-r));
exp_ra1=(1-exp_ra)/(2*palpha/(1-r));

hk=abs((va_n+kesy1)/(va+kesy1))^(-(1+r)/(1-r));
gk=2/(pbeta+sqt_baw)*(1-exp_ra*hk);
gk1=2/(pbeta+sqt_baw)*(dt-exp_ra1*hk);

%计算矩阵
vx=zeros(3,3);%方向分量
vx(1,2)=-v(3); vx(1,3)=v(2);
vx(2,1)=v(3);  vx(2,3)=-v(1);
vx(3,1)=-v(2); vx(3,2)=v(1);

A(1:3,4:6)=eye(3)*(exp_ra1*hk);
A(1:3,7)=gk1*v/va;
A(1:3,8:10)=-(exp_ra1*hk+gk1*w/va)*dt*vx;
A(4:6,4:6)=eye(3)*(exp_ra*hk);
A(4:6,7)=gk*v/va;
A(4:6,8:10)=-(exp_ra*hk+gk*w/va)*dt*vx;

xn=A*x;
end

function xn=discreteModel2(x,palpha,pbeta,dt)
%1.数据准备
n=length(x);
A=eye(n,n);

%速度
v=x(4:6);   %速度分量
va=sqrt(v'*v);
if va<0.00001
    va=0.00001;
end

%功率
w=x(7);

%2.计算下一时刻速度
ra=w/pbeta;

va_n=va;
f2=va+ra*log(abs(va-ra))-pbeta*dt;  %前一时刻系数
for iter=1:10
    f1=va_n+ra*log(abs(va_n-ra));
    f=f1-f2;
    df=va_n/(va_n-ra);
    va_n=va_n-(va_n/ra)*f/df;
end

%3.状态转移矩阵
%基本参数
exp_rb=exp(-pbeta*dt/ra);%指数量
exp_rb1=(1-exp_rb)*ra/pbeta;%指数量一阶积分

hk=exp((va-va_n)/ra);
gk=(1-exp_rb*hk)/(pbeta*va);
gk1=(dt-exp_rb1*hk)/(pbeta*va);

%矩阵
vx=zeros(3,3);%方向分量
vx(1,2)=-v(3); vx(1,3)=v(2);
vx(2,1)=v(3);  vx(2,3)=-v(1);
vx(3,1)=-v(2); vx(3,2)=v(1);

A(1:3,4:6)=(exp_rb1*hk)*eye(3);
A(1:3,7)=gk1*v;
A(1:3,8:10)=-(exp_rb1*hk+gk1*w)*dt*vx;

A(4:6,4:6)=(exp_rb*hk)*eye(3);
A(4:6,7)=gk*v;
A(4:6,8:10)=-(exp_rb*hk+gk*w)*dt*vx;

xn=A*x;
end

function xn=DPTmodel2(x,palpha,pbeta,dt)
%数据准备
n=length(x);
Fy=eye(n,n);

%速度
v=x(4:6);   %速度分量
va=sqrt(v'*v);
if va<0.001
    va=0.001;
end

%功率
w=x(7);

%计算速度参数
exp_at=exp(-palpha*dt);
exp_at1=(1-exp_at)/palpha;
exp_at2=(dt-exp_at1)/palpha;

gk=exp_at-pbeta/va*exp_at1;
gk1=exp_at1-pbeta/va*exp_at2;
hk=exp_at1/(va^2+dt*w/4);
hk1=exp_at2/(va^2+dt*w/4);
fiyk=dt*(gk+hk*w);
fiyk1=dt*(gk1+hk1*w);

%计算矩阵
vx=zeros(3,3);%方向分量
vx(1,2)=-v(3); vx(1,3)=v(2);
vx(2,1)=v(3);  vx(2,3)=-v(1);
vx(3,1)=-v(2); vx(3,2)=v(1);

Fy(1:3,4:6)=eye(3)*gk1;
Fy(1:3,7)=hk1*v;
Fy(1:3,8:10)=-fiyk1*vx;
Fy(4:6,4:6)=eye(3)*gk;
Fy(4:6,7)=hk*v;
Fy(4:6,8:10)=-fiyk*vx;

xn=Fy*x;
end
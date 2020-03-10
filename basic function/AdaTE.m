%基于 adaptive trajectory estimation 的轨迹估计函数
%使用DPT模型
function [As_up,Qs_up,Rzs_up,cut_t,Xs_up,dX_up,preX]=AdaTE(As,Qs,Rzs_ob,Rzs,H,cut_t,Xs,dX,Zs,sitar,time)
%@As 状态转移矩阵
%@Qs 观测协方差
%@Rzs_ob 观测方差的初始值
%@Rzs 修正后的观测方差
%@H 观测矩阵
%@cut_t 截断时刻
%@Xs 轨迹状态
%@dX 转移误差
%@Zs 观测
%@sitar 参数
%@time 时间

% 参数准备
palpha=sitar.alpha;
pbeta=sitar.beta;
Da=sitar.Da;
Dt=sitar.Dt;

[~,nx]=size(Xs);
[~,nz]=size(Zs);

if(nz<=nx)
    preX=-1;    %观测信息不足以滤波，报错误
    disp('error:观测信息不足，无法完成滤波，请检查观测信息');
    return;
end

% 1.更新状态转移矩阵
As_up=[As(1:nx);cell(1)];
Qs_up=[Qs(1:nx);cell(1)];
for i=cut_t+1:nx-1
    dt=time(i+1)-time(i);
    As_up{i}=MotionModal2(Xs(:,i),dt,palpha,pbeta);
end

% 2.更新最后一个协方差
if(nx>0)
    dt=time(nx+1)-time(nx);
    [As_up{nx},Qs_up{nx}]=transModel(Xs(:,nx),dt,palpha,pbeta,Da,Dt);
end

% 3.轨迹预测
if(nx==0)
    preX=zeros(10,1);
    preX(1:3,1)=Zs(1:3,1);
    preX(4:6,1)=(rand(3,1)-0.5);
else
    preX=As_up{nx}*Xs(:,nx);
end

% 4.更新观测方差
Rzs_up=observCovUpdate(Rzs_ob,Rzs,H,cut_t,[Xs,preX],Zs);

% 5.计算新的轨迹
Xs_up=MAPestimation(As_up,Qs_up,Rzs_up,H,cut_t,Xs,Zs,time);

% 6.更新状态转移方差
dX_up=stransError(As_up,Xs_up,dX,cut_t);  %更新转移误差
Qs_up=transCovUpdate(Qs_up,dX_up,time,cut_t);

% 7.截断时刻更新
Fiy=diag([1,1,1,1/4,1/4,1/4,1/9,1/9,1/9,1/9]);
cut_t=cutTimeUpdate([Xs,preX],Xs_up,Fiy,cut_t);

end

function A=MotionModal(x,dt,palpha,pbeta)
%  计算状态转移矩阵A,偏转矩阵F
%  输入： 上一时刻状态向量x 时间间隔dt 阻尼因子 palpha,pbeta
%  输出：状态转移矩阵A

n=length(x);
A=eye(n,n);

w=x(7);
if w<0
    w=0;
end

%基本参数
v=x(4:6);   %速度分量
v2=v'*v;
if(v2<0.000001)
    v2=0.000001;
end
va=sqrt(v2);

% 1.计算状态转移与偏移矩阵
if w<pbeta^2/(palpha*4)
    r=pbeta/sqrt(pbeta^2+4*palpha*w);
    r=2*palpha/(1+r);
    
    exp_r=exp(-r*dt);%指数量
    exp_r1=(1-exp_r)/r;%指数量一阶积分
    texp_r=dt*exp_r;%时间指数量
    texp_r1=(exp_r1-texp_r)/r;%时间指数量
    
    peta_v=pbeta/va;
    
    %位移
    for i=1:3
        A(i,i+3)=exp_r1-dt^2*peta_v/2;
    end
    A(1:3,8:10)=dt^3*peta_v/3-texp_r1;
    
    % 速度增量
    for i=4:6
        A(i,i)=exp_r-dt*peta_v;
    end
    A(4:6,8:10)=dt^2*peta_v-texp_r;
else
    exp_a=exp(-palpha*dt);    %指数量
    exp_a1=(1-exp_a)/palpha;%指数量一阶积分
    texp_a=dt*exp_a;    %时间指数量
    texp_a1=(exp_a1-texp_a)/palpha;%时间指数一阶积分
    
    para_gt=v2+dt*w/2;
    para_w=v/para_gt;       %功率分量
    para_ct=zeros(3,3);%方向分量
    para_ct(1,2)=-v(3); para_ct(1,3)=v(2);
    para_ct(2,1)=v(3);  para_ct(2,3)=-v(1);
    para_ct(3,1)=-v(2); para_ct(3,2)=v(1);
    
    para_ct=para_ct*(1+dt*w/para_gt);
    
    for i=1:3
        A(i,i+3)=exp_a1;
    end
    
    A(1:3,7)=texp_a1*para_w;
    A(1:3,8:10)=-texp_a1*para_ct;
    
    % 速度增量
    for i=4:6
        A(i,i)=exp_a;
    end
    
    A(4:6,7)=texp_a*para_w;
    A(4:6,8:10)=-texp_a*para_ct;
end
end

function A=MotionModal2(x,dt,palpha,pbeta)
%  计算状态转移矩阵A,使用偏导数近似计算
%  输入： 状态向量x,xp 下一个速度vn 时间间隔dt 阻尼因子 palpha,pbeta
%  输出：状态转移矩阵A

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
if w<-pbeta^2/(4*palpha);
    w=0.000001-pbeta^2/(4*palpha);
end

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

if isnan(A(1,6))
    A(1,6);
end

end

function A=MotionModal3(x,dt,palpha,pbeta)
%  计算状态转移矩阵A,使用偏导数近似计算
%  输入： 状态向量x,xp 下一个速度vn 时间间隔dt 阻尼因子 palpha,pbeta
%  输出：状态转移矩阵A

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
if abs(w)<0.01
    w=0.01*w/abs(w);
end

%2.计算下一时刻速度
ra=w/pbeta;

va_n=va;
f2=va+ra*log(abs(va-ra))-pbeta*dt;  %前一时刻系数
for iter=1:10
    f1=va_n+ra*log(abs(va_n-ra));
    f=f1-f2;
    df=va_n/(va_n-ra);
    
    rt=abs(va_n/ra);
    if rt>1
        rt=1/rt;
    end
    va_n=va_n-rt*f/df;
    va_n=abs(va_n);
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

if isnan(A(1,6))
    A(1,6);
end

end


function [A,Q]=transModel(x,dt,palpha,pbeta,Da,Dt)
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
if w<-pbeta^2/(4*palpha);
    w=0.000001-pbeta^2/(4*palpha);
end

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

% 4.协方差
Q=zeros(n,n);

%叉乘矩阵
vxDt=vx*Dt;
vxDtvxt=vxDt*Dt';
Vvt=v*v';

Q(1:3,1:3)=Da*dt*(gk1/va)^2/3*Vvt+dt^3*(exp_ra1*hk+gk1*w/va)^2/5*vxDtvxt;
Q(1:3,4:6)=Da*dt*(gk1*gk/va^2)/2*Vvt+dt^3*(exp_ra1*hk+gk1*w/va)*(exp_ra*hk+gk*w/va)/4*vxDtvxt;
Q(1:3,7)=Da*dt*gk1/(va*2)*v;
Q(1:3,8:10)=-dt*(exp_ra1*hk+gk1*w/va)/2*vxDt;

Q(4:6,4:6)=Da*dt*(gk/va)^2*Vvt+dt^3*(exp_ra*hk+gk*w/va)^2/3*vxDtvxt;
Q(4:6,7)=Da*dt*gk/(va)*v;
Q(4:6,8:10)=-dt*(exp_ra*hk+gk*w/va)*vxDt;

Q(7,7)=Da;
Q(8:10,8:10)=Dt;

%对称化
for i=1:n
    for j=i+1:n
        Q(j,i)=Q(i,j);
    end
end

for i=1:n
    Q(i,i)=Q(i,i)*n/2;
end

end

function [A,Q]=transModel2(x,dt,palpha,pbeta,Da,Dt)
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
    
    rt=abs(va_n/ra);
    if rt>1
        rt=1/rt;
    end
    va_n=va_n-rt*f/df;
    va_n=abs(va_n);
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

% 4.协方差
Q=zeros(n,n);

%叉乘矩阵
vxDt=vx*Dt;
vxDtvxt=vxDt*Dt';
Vvt=v*v';

fiy=(exp_rb*hk+gk*w)*dt;

Q(1:3,1:3)=Da*dt*(gk1)^2/3*Vvt+dt*(fiy*dt)^2/5*vxDtvxt;
Q(1:3,4:6)=Da*dt*gk1*gk^2/2*Vvt+dt*(fiy^2*dt)/4*vxDtvxt;
Q(1:3,7)=Da*dt*gk1*v/2;
Q(1:3,8:10)=-dt*(fiy*dt)/3*vxDt;

Q(4:6,4:6)=Da*dt*(gk)^2*Vvt+dt*(fiy)^2/3*vxDtvxt;
Q(4:6,7)=Da*dt*gk*v;
Q(4:6,8:10)=-dt*(fiy)/2*vxDt;

Q(7,7)=Da;
Q(8:10,8:10)=Dt;

%对称化
for i=1:n
    for j=i+1:n
        Q(j,i)=Q(i,j);
    end
end

for i=1:n
    Q(i,i)=Q(i,i)*n/2;
end

end

% 4.更新观测方差
function Rzs=observCovUpdate(Rzs_ob,Rzs,H,cut_t,Xs,Zs)
[nr,lz]=size(Rzs);
mz=size(Zs,1);

if(nr<1)
    return;
end

% 计算观测误差
dz=zeros(mz,nr-cut_t,lz);
rz=zeros(nr-cut_t,lz);
for i=cut_t+1:nr
    id=i-cut_t;
    for j=1:lz
        cur_dz=H*Xs(:,i)-Zs(:,i,j);
        dz(:,id,j)=cur_dz;
        rz(id,j)=cur_dz'/Rzs_ob{i,j}*cur_dz;
    end
end

h_rz=min([3*sqrt(mean(mean(rz.^2))),16]);

for i=cut_t+1:nr
    id=i-cut_t;
    for j=1:lz
        if(rz(id,j)>h_rz)
            Rzs{i,j}=Rzs_ob{i,j}/2+dz(:,id,j)*dz(:,id,j)';
        else
            Rzs{i,j}=Rzs_ob{i,j};
        end
    end
end

end

% 5.计算新的轨迹
function Xs_up=MAPestimation(As,Qs,Rzs,H,cut_t,Xs,Zs,time)
%  trajectory MAP from observations
%  输入：离散状态转移矩阵集合As,状态转移协方差矩阵集合Qs,
%       观测协方差矩阵Rzs,观测矩阵H,截断时刻cut_t,上一时刻轨迹估计Xs,观测Zs
%  输出：状态轨迹X_re

n=length(As);  %状态时刻数
[nz,lz]=size(Rzs);%观测时刻数及同一时刻观测个数
if(nz<n)
    disp('error->in maxLikelihoodFilter 观测个数过少，无法求解');
    return;
end
if(n==0)
    disp('error->in maxLikelihoodFilter 状态转移矩阵数目小于1 无法求解');
    return;
end

%只有少许状态，直接计算初始值
if n<=5
    Xs_up=0.1*(rand(10,n)-0.5);
    for i=1:n
        Xs_up(1:3,i)=positionFromObserv(Zs(:,i,:),Rzs(i,:));
        Xs_up(7,i)=abs(Xs_up(7,i))*10;
        if i>=2
            Xs_up(4:6,i-1)=Xs_up(4:6,i-1)+(Xs_up(1:3,i)-Xs_up(1:3,i-1))/(time(i)-time(i-1));%速度
        end
    end
    
    if n>1
        Xs_up(4:6,n)=Xs_up(4:6,n)+Xs_up(4:6,n-1);
    end
    return;
end

if n>80
    cut_t=max(cut_t,2);
end

%1.包含多个状态的迭代求解
%(1)计算方差的逆
n_use=n-cut_t;

Qs_inv=cell(n_use,1);
Rs_inv=cell(n_use,lz);

for i=1:n_use
    t=cut_t+i;  %时间戳
    
    Q=Qs{t};    %状态转移协方差
    Qs_inv{i}=inv(Q);   %逆矩阵
    
    for j=1:lz
        Rs_inv{i,j}=inv(Rzs{t,j});  %观测协方差逆
    end
end

if(cut_t>0) %前一个协方差阵求逆
    Qp_inv=inv(Qs{cut_t});
end

% (2)系数矩阵与向量的各个分量
m=size(Xs,1);

Mc=cell(n_use,1);   %M的对角元
Ms=cell(n_use,1);   %M的非对角元
b=zeros(m,n_use);   %每列同一个时刻，行方向为时间轴

% 第一个时刻
t=cut_t+1;
xm_Mat=As{t}'*Qs_inv{1};
Ms{1}=xm_Mat;
b_cr=Rs_inv{1,1}*Zs(:,t,1);
sum_Rs=Rs_inv{1,1};
for j=2:lz
    sum_Rs=sum_Rs+Rs_inv{1,j};
    b_cr=b_cr+Rs_inv{1,j}*Zs(:,t,j);
end
Mc{1}=H'*sum_Rs*H+xm_Mat*As{t};
b(:,1)=H'*b_cr;

if(cut_t>0)    %截断时刻生效
    Mc{1}=Mc{1}+Qp_inv; %t==cut_t
    b(:,1)=b(:,1)+Qp_inv*As{cut_t}*Xs(:,cut_t);
end

% 后续时刻
for i=2:n_use-1
    t=cut_t+i;  %时间戳
    xm_Mat=As{t}'*Qs_inv{i};
    Ms{i}=xm_Mat;
    sum_Rs=Rs_inv{i,1};
    b_cr=Rs_inv{i,1}*Zs(:,t,1);
    
    for j=2:lz
        sum_Rs=sum_Rs+Rs_inv{i,j};
        b_cr=b_cr+Rs_inv{i,j}*Zs(:,t,j);
    end
    
    Mc{i}=H'*sum_Rs*H+xm_Mat*As{t}+Qs_inv{i-1}; %n_use<t<k
    b(:,i)=H'*b_cr;
end

%t==k
t=cut_t+n_use;
sum_Rs=Rs_inv{n_use,1};
b_cr=Rs_inv{n_use,1}*Zs(:,t,1);

for j=2:lz
    sum_Rs=sum_Rs+Rs_inv{n_use,j};
    b_cr=b_cr+Rs_inv{n_use,j}*Zs(:,t,j);
end
Mc{n_use}=H'*sum_Rs*H+Qs_inv{n_use-1};
b(:,n_use)=H'*b_cr;

%2.根据观测值校正轨迹
% (1)正序消元
for i=1:n_use-1
    cg=Ms{i}'/Mc{i};
    Mc{i+1}=Mc{i+1}-cg*Ms{i};
    b(:,i+1)=b(:,i+1)+cg*b(:,i);
end

% (3)逆序推算
X_up=zeros(m,n_use);	%X更新部分
X_up(:,n_use)=Mc{n_use}\b(:,n_use);

if isnan(X_up(1,n_use))
    t=cut_t+n_use;
    X_up(:,n_use)=Xs(:,t-1);
    X_up(1:3,n_use)=positionFromObserv(Zs(:,t,:),Rzs(t,:));
    X_up(4:6,n_use)=(X_up(1:3,n_use)-X_up(1:3,n_use-1))/(time(t)-time(t-1));
end

for i=n_use-1:-1:1
    X_up(:,i)=Mc{i}\(b(:,i)+Ms{i}*X_up(:,i+1));
    if isnan(X_up(1,i))
        X_up(:,1:i)=Xs(:,cut_t+1:cut_t+i)+0.04*(rand(10,i)-0.5);
        break;
    end
end

% 输出更新
Xs_up=[Xs(:,1:cut_t),X_up];
end

% 计算状态转移误差
function dX_re=stransError(As,Xs,dX,cut_t)
n=length(As);
% 状态转移误差
dX_re=[dX,zeros(10,1)];
for i=cut_t+1:n-1
    dX_re(:,i)=As{i}*Xs(:,i)-Xs(:,i+1);
end
end

% 状态转移方差更新函数
function Qs=transCovUpdate(Qs,delta_x,time,cut_t)
% 输入：状态转移矩阵Qs 当前轨迹状态转移误差delta_x 采样时刻time 
%       截断时刻cut_t 当前机动tm

nt=length(Qs);
if(nt<2)
    return;
end

time_use=time(nt)-time(cut_t+1);
dt=time(nt)-time(nt-1);

for i=cut_t+1:nt-1
    dd=delta_x(:,i);%-delta_pr(:,i);
    Qs{i}=(time_use*Qs{i}+(dt*dd)*dd')/(time_use+dt);
end

end

% 截断时刻更新
function cut_t=cutTimeUpdate(Xs_pre,Xs_ev,Fiy,cut_tp)
%  输入：状态轨迹预测Xs_pre 状态轨迹估计Xs_ev 权重矩阵Fiy 前一个截断时刻cut_tp
%  输出：截断时刻cut_t

[~,n]=size(Xs_pre);
ang=0.000001;        %误差限

sum_eer=0;

for i=cut_tp+1:2:n-40   %隔两个时刻一次采样
    error=Xs_pre(:,i)-Xs_ev(:,i);
    if sum_eer==0
        sum_eer=error'*Fiy*error;
    else
        sum_eer=0.95*sum_eer+0.05*error'*Fiy*error;
    end
    if(sum_eer>ang)    %平均误差大于阈值
        cut_t=max(i-20,0);
        return;
    end
end
cut_t=max(n-50,0);
end

function p=positionFromObserv(zs,Rz)
lz=length(Rz);%观测时刻数及同一时刻观测个数

p=zeros(3,1);
sum_w=zeros(3,3);

for j=1:lz
    iv_Rz=inv(Rz{j});
    p=p+iv_Rz*zs(:,j);
    sum_w=sum_w+iv_Rz;
end

p=sum_w\p;
end

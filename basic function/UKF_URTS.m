%UKF-URTS for trajectory estimation
%
%input:
% observation:obs
% observe matrix:H,
% covariance of observation: Rzs_ob
% state covariance at biginning: Qs,
% transition covariance: Qt,
% previous estimation: Xs,
% real time localization: Xs_real,
% some other parameters: palpha,pbeta,time,cut_t,dynamic
%
%output:
% updated trajectory: Xs_up
% updated covariances: Qs_up

function [Xs,Qs]=UKF_URTS(Zs,H,Rzs,Qs,Qt,Xs,Xs_real,palpha,pbeta,time,cut_t,dynamic)

[mx,nx]=size(Xs_real);
nz=size(Zs,2);
cur_t=nx+1;

if(nz<=nx)
    disp('error:观测信息不足，无法完成滤波，请检查观测信息');
    return;
end

% 若存在多个观测，校正观测
lz=size(Rzs,2);

if(lz>1)    %观测数量大于1个，取加权均值
    sinv_Rzs=zeros(size(Rzs{cur_t,1}));
    cur_z=zeros(length(Rzs{cur_t,1}),1);
    
    for j=1:lz
        inv_R=inv(Rzs{cur_t,j});
        sinv_Rzs=sinv_Rzs+inv_R;
        cur_z=cur_z+inv_R*Zs(:,cur_t,j);
    end
    
    Rz=inv(sinv_Rzs);
    cur_z=sinv_Rzs\cur_z;
else
    Rz=Rzs{cur_t,1};
    cur_z=Zs(:,cur_t,1);
end

% 1.UKF滤波更新最后一个状态
Xs=[Xs,zeros(mx,1)];
Qs=[Qs;cell(1)];
[Xs(:,nx+1),Qs{nx+1}]=UnscentedKalmanFilter(cur_z,H,Rz,Qt,Xs(:,nx),Qs{nx},palpha,pbeta,time(nx+1)-time(nx),dynamic);

% 2.URTS平滑更新之前状态
Xs=URTSsmoother(Xs,Xs_real,Qt,Qs,palpha,pbeta,time,cut_t,dynamic);
end

%% UKF滤波器
function [x_ev,Q_ev]=UnscentedKalmanFilter(z,H,Rz,Qt,xk,Q_pre,palpha,pbeta,dt,choose)

% 基本参数
nx=length(xk);
a=1;
beita=3;
kap=5;
lambda=a^2*(nx+kap)-nx;
c=nx+lambda;

wm=ones(nx*2+1,1)/(2*c);  %权重
wm(nx*2+1)=lambda/c;    %均值权重
wc=wm;  %权重
wc(nx*2+1)=lambda/c+1-a^2+beita; %均值权重

% 1.k时刻Sigma采样
xk_sigma=zeros(nx,2*nx+1);
root_Q=sqrtm(c*Q_pre);

for i=1:nx
    xk_sigma(:,i)=xk+root_Q(:,i);
    xk_sigma(:,nx+i)=xk-root_Q(:,i);
end
xk_sigma(:,nx*2+1)=xk;

% 2.采样状态转移预测
x_predict=zeros(nx,2*nx+1);
if(choose=='Acc')   %加速度模型
    for i=1:2*nx+1
        x_predict(:,i)=accModel(xk_sigma(:,i),pbeta,dt);
    end
else            %功率模型
    for i=1:2*nx+1
        x_predict(:,i)=DPTmodel2(xk_sigma(:,i),palpha,pbeta,dt);
    end
end

% 3.计算预测转移均值与方差
mx_pre=zeros(nx,1);
Px_pre=Qt;

for i=1:2*nx+1
    mx_pre=mx_pre+x_predict(:,i)*wm(i);
end

for i=1:2*nx+1
    dx=x_predict(:,i)-mx_pre;
    Px_pre=Px_pre+(wc(i)*dx)*dx';
end

% 4.预测值sigma采样
root_Ppre=sqrtm(c*Px_pre);

for i=1:nx
    xk_sigma(:,i)=mx_pre+root_Ppre(:,i);
    xk_sigma(:,nx+i)=mx_pre-root_Ppre(:,i);
end
xk_sigma(:,nx*2+1)=mx_pre;

% 5.采样观测预测
nz=length(z);
z_pre=zeros(nz,nx*2+1);
for i=1:nx*2+1
    z_pre(:,i)=H*xk_sigma(:,i);
end

mz_pre=zeros(nz,1);
Sz=Rz;  %观测方差加上转移方差
Cxz=zeros(nx,nz);
for i=1:nx*2+1
    mz_pre=mz_pre+z_pre(:,i)*wm(i);
end

for i=1:nx*2+1
    dx=x_predict(:,i)-mx_pre;
    dz=z_pre(:,i)-mz_pre;
    
    Sz=Sz+(wc(i)*dz)*dz';
    Cxz=Cxz+(wc(i)*dx)*dz';
end

% 6.状态校正
Ka=Cxz/Sz;
x_ev=mx_pre+Ka*(z-mz_pre);
Q_ev=Px_pre-Ka*Cxz';    %Ka*Sz*Ka'=Ka*Sz*(Cxz/Sz)'=Ka*Cxz'
end

%% URTS平滑器
function Xs=URTSsmoother(Xs,Xs_real,Qt,Qs,palpha,pbeta,time,cut_t,choose)
[L,nx]=size(Xs);

% 基本参数
a=1;
beita=3;
kap=5;
lambda=a^2*(L+kap)-L;
c=L+lambda;

wm=ones(L*2+1,1)/(2*c);  %权重
wm(L*2+1)=lambda/c;    %均值权重
wc=wm;  %权重
wc(L*2+1)=lambda/c+1-a^2+beita; %均值权重

% 滤波
xk_sigma=zeros(L,2*L+1);
x_predict=zeros(L,2*L+1);

for i=nx-1:-1:cut_t
    dt=time(i+1)-time(i);
    
    % 1.sigma采样
    root_Q=sqrtm(c*Qs{i});
    xk=Xs_real(:,i);
    
    for j=1:L
        xk_sigma(:,j)=xk+root_Q(:,j);
        xk_sigma(:,L+j)=xk-root_Q(:,j);
    end
    xk_sigma(:,L*2+1)=xk;
    
    % 2.采样状态转移预测
    if(choose=='Acc')   %加速度模型
        for j=1:2*L+1
            x_predict(:,j)=accModel(xk_sigma(:,j),pbeta,dt);
        end
    else            %功率模型
        for j=1:2*L+1
            x_predict(:,j)=DPTmodel2(xk_sigma(:,j),palpha,pbeta,dt);
        end
    end
    
    % 3.计算预测转移均值与方差
    mx_pre=zeros(L,1);  %预测均值
    Px_pre=Qt;          %预测方差
    Cx_pre=zeros(L,L);  %预测协方差
    
    for j=1:2*L+1
        mx_pre=mx_pre+x_predict(:,j)*wm(j);
    end

    for j=1:2*L+1
        dx=x_predict(:,j)-mx_pre;
        dx_k=xk_sigma(:,j)-xk;
        
        Px_pre=Px_pre+(wc(j)*dx)*dx';
        Cx_pre=Cx_pre+(wc(j)*dx_k)*dx';
    end
    
    % 4.平滑校正
    Xs(:,i)=xk+Cx_pre*(Px_pre\(Xs(:,i+1)-mx_pre));
end
end

%% 运动模型
% 加速度模型
function x_new=accModel(x,pbeta,dt)
x_new=x;

% 速度->位置
dt_2=dt^2/2;
x_new(1:3)=x_new(1:3)+(dt-pbeta*dt_2)*x(4:6)+dt_2*x(7:9);

% 加速度->速度 及阻尼
x_new(4:6)=(1-pbeta*dt)*x(4:6)+dt*x(7:9);
end

% 功率模型
function x_new=DPTmodel(x,palpha,pbeta,dt)
%数据准备
n=length(x);
Fy=eye(n,n);

%速度
v=x(4:6);   %速度分量
va=sqrt(v'*v);
if va<0.00001
    va=0.00001;
end

%功率
w=x(7);
if abs(w)<0.001
    w=w/abs(w)*0.001;
end
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

Fy(1:3,4:6)=eye(3)*(exp_ra1*hk);
Fy(1:3,7)=gk1*v/va;
Fy(1:3,8:10)=-(exp_ra1*hk+gk1*w/va)*dt*vx;
Fy(4:6,4:6)=eye(3)*(exp_ra*hk);
Fy(4:6,7)=gk*v/va;
Fy(4:6,8:10)=-(exp_ra*hk+gk*w/va)*dt*vx;

x_new=Fy*x;
end

function x_new=DPTmodel2(x,palpha,pbeta,dt)
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

%计算矩阵
ct=x(8:10);
cx=zeros(3,3);%方向分量
cx(1,2)=-ct(3); cx(1,3)=ct(2);
cx(2,1)=ct(3);  cx(2,3)=-ct(1);
cx(3,1)=-ct(2); cx(3,2)=ct(1);

vn=(exp_at+exp_at1*(w/va^2-pbeta/va))*expm(dt*cx)*v;
pn=x(1:3)+dt*(v+vn)/2;

x_new=[pn;vn;x(7:10)];
end

%EKF-RTS for trajectory estimation
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

function [Xs,Qs]=EKF_RTS(Zs,H,Rzs,Qs,Qt,Xs,Xs_real,palpha,pbeta,time,cut_t,choose)

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

% 1.kalman滤波更新最后一个状态
Xs=[Xs,zeros(mx,1)];
Qs=[Qs;cell(1)];
[Xs(:,cur_t),Qs{cur_t}]=ExtendedKalmanFilter(cur_z,H,Rz,Qt,Xs_real(:,nx),Qs{nx},palpha,pbeta,time(cur_t)-time(nx),choose);

% 2.RTS平滑更新之前状态
As=cell(nx,1);
if(choose=='Acc')	%加速度模型
    for i=cut_t:nx
        dt=time(i+1)-time(i);
        As{i}=accModel(Xs(:,i),pbeta,dt);
    end
else
    for i=cut_t:nx
        dt=time(i+1)-time(i);
        As{i}=DPTmodel2(Xs(:,i),palpha,pbeta,dt);
    end
end

Xs=RTSsmoother(Xs,Xs_real,As,Qt,Qs,cut_t);

end

%% EKF滤波器
function [x_ev,Q_ev]=ExtendedKalmanFilter(z,H,Rz,Qt,x_pre,Q_pre,palpha,pbeta,dt,choose)
% 状态转移
if(choose=='Acc')   %加速度模型
    Fy=accModel(x_pre,pbeta,dt);
else
    Fy=DPTmodel2(x_pre,palpha,pbeta,dt);
end

% 预测
x_ev=Fy*x_pre;
Q_ev=Fy*Q_pre*Fy'+Qt;

% 校正
K=Q_ev*H'/(H*Q_ev*H'+Rz);
x_ev=x_ev+K*(z-H*x_ev);
Q_ev=Q_ev-K*(H*Q_ev);
end

%% RTS平滑器
function Xs=RTSsmoother(Xs,Xs_real,As,Qt,Qs,cut_t)
%  输入：状态轨迹Xs 实时估计Xs_real 状态漂移mrs 状态转移As 转移方差Qt 状态方差Qs 截断时刻cut_t
[~,nx]=size(Xs);

for i=nx-1:-1:cut_t
    % 预测值
    A=As{i};
    x_pre=A*Xs_real(:,i);
    P=Qs{i};
    Q_pre=A*P*A'+Qt;
    
    % 后验校正
    K=P*A'/Q_pre;
    Xs(:,i)=Xs_real(:,i)+K*(Xs(:,i+1)-x_pre);
end

end

%% 运动模型
% Acc model
function Fy=accModel(x,pbeta,dt)
n=length(x);
Fy=eye(n);

% 速度->位置
for i=1:3
    Fy(i,i+3)=dt;
    Fy(i,i+6)=dt^2/2;
end

% 加速度->速度 及阻尼
for i=4:6
    Fy(i,i+3)=dt;
    Fy(i,i)=Fy(i,i)-pbeta*dt;
end
end

%DPT model
function Fy=DPTmodel(x,palpha,pbeta,dt)
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

end

function Fy=DPTmodel2(x,palpha,pbeta,dt)
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

end
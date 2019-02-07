%基于 online adaptive optimization 的轨迹估计函数
%使用DPT模型
%用于路径规划
function [As_up,Qs_up,cut_t,Xs_up,dX_up,preX]=OAOestimation_PathPlan(As,Qs,Rzs_kp,Rzs_br,H,cut_t,Xs,dX,Zkp,Zbr,sitar,time,dt)
%@As 状态转移矩阵
%@Qs 观测协方差
%@Rzs_kp 关键点的允许偏离方差
%@Rzs_br 障碍物的尺寸对应方差
%@H 观测矩阵
%@cut_t 截断时刻(当前时刻)
%@Xs 轨迹状态
%@dX 转移误差
%@Zkp 路径关键点坐标
%@Zbr 路径障碍物坐标
%@sitar 参数
%@time 时间(到达各个key point 时间)
%@dt 轨迹采样时间间隔

% 参数准备
palpha=sitar.alpha;
pbeta=sitar.beta;
Da=sitar.Da;
Dt=sitar.Dt;

[~,nx]=size(Xs);
[~,nzk]=size(Zkp);  %key point个数

if(nzk<=1)
    preX=-1;    %观测信息不足以滤波，报错误
    disp('error:观测信息不足，无法完成滤波，请检查观测信息');
    return;
end

%选取后续2个关键点
Zxp=zeros(3,2);
nxp=zeros(2,1);
nxt=zeros(2,1);
fut_p=time>cut_t;
nn=0;
for i=1:nzk
    if fut_p(i)
        nn=nn+1;
        Zxp(:,nn)=Zkp(:,i);
        nxp(nn)=i;
        nxt(nn)=time(i);
        
        if nn==2
            break;
        end
    end
end

%状态初始化
lx=size(Xs(:,nxt(1)),1);
if sum(Xs(:,nxt(1))==0)    %初始化条件
    cur_x=zeros(lx,1);
    cur_x(1:3)=Zxp(:,1);
    cur_x(4:6)=(Zxp(:,2)-Zxp(:,1))/(nxt(2)-nxt(1));
    v2=cur_x(4:6)'*cur_x(4:6);
    va=sqrt(v2)';
%     cur_x(7)=sqrt(v2)*pbeta+v2*palpha;
%     cur_x(8:10)=cur_x(4:6)/v2;
    cur_x(7:9)=cur_x(4:6)*(pbeta/va+va*palpha);
    Xs(:,nxt(1))=cur_x+[0;0;0;rand(lx-3,1)*0.1-0.05];
end

if sum(Xs(:,nxt(2))==0)    %初始化条件
    cur_x=Xs(:,nxt(1));
    cur_x(1:3)=Zxp(:,2);
    Xs(:,nxt(2))=cur_x+[0;0;0;rand(lx-3,1)*0.1-0.05];
end

% 1.更新状态转移矩阵
As_up=[As(1:nx);cell(nxt(2)-nx)];
Qs_up=[Qs(1:nx);cell(nxt(2)-nx)];
for i=cut_t+1:nxt(2)-1
    As_up{i}=AccModel(Xs(:,i),dt,palpha,pbeta);
%     As_up{i}=MotionModal2(Xs(:,i),dt,palpha,pbeta);
    preX=As_up{i}*Xs(:,i);  %预测
    
    if sum(Xs(:,i+1)==0)
        Xs(:,i+1)=preX+[0;0;0;rand(lx-3,1)*0.01-0.005];
    end
end

% 2.更新最后一个协方差
if(cut_t==0||cut_t==time(nxp(1)-1))    %刚到达新的点
    for i=cut_t+1:nxt(2)-1
%         [As_up{i},Qs_up{i}]=transModel(Xs(:,i),dt,palpha,pbeta,Da,Dt);
        [As_up{i},Qs_up{i}]=transModelAcc(Xs(:,i),dt,palpha,pbeta,Da);
        
        preX=As_up{i}*Xs(:,i);  %预测
        
        if sum(Xs(:,i+1)==0)
            Xs(:,i+1)=preX+[0;0;0;rand(lx-3,1)*0.01-0.005];
        end
    end
end

% 4.更新观测方差
%Rzs_up=observCovUpdate(Rzs_ob,Rzs,H,cut_t,[Xs,preX],Zs);

% 5.计算新的轨迹

if cut_t==0
    Xs_up=Xs;%初始化，约束不足无法求解，直接退出
else
    Xs_up=MAPestimationPP(As_up,Qs_up,H,cut_t,Xs,Rzs_kp(nxp(1):nxp(2)),Rzs_br,Zxp,Zbr,nxt);
end

% 6.更新状态转移方差
dX_up=stransError(As_up,Xs_up,dX,cut_t);  %更新转移误差
Qs_up=transCovUpdate(Qs_up,dX_up,cut_t,dt);

% 7.截断时刻更新
% Fiy=diag([1,1,1,1/4,1/4,1/4,1/9,1/9,1/9,1/9]);
% cut_t=cutTimeUpdate([Xs,preX],Xs_up,Fiy,cut_t);
end

function A=AccModel(x,dt,palpha,pbeta)
%  计算状态转移矩阵A
%  输入： 上一时刻状态向量x 时间间隔dt 阻尼因子 palpha,pbeta
%  输出：状态转移矩阵A

%数据准备
n=length(x);
A=eye(n,n);

%基本参数
v=x(4:6);   %速度分量
v2=v'*v;
if(v2<0.0001)
    v2=0.0001;
end
va=sqrt(v2);

exp_at=exp(-palpha*dt);
exp_at1=(1-exp_at)/palpha;
exp_at2=(dt-exp_at1)/palpha;

A(1:3,4:6)=eye(3)*(exp_at1-exp_at2*pbeta/va);
A(1:3,7:9)=eye(3)*exp_at2/va;
A(4:6,4:6)=eye(3)*(exp_at-exp_at1*pbeta/va);
A(4:6,7:9)=eye(3)*exp_at1/va;

% A=eye(n,n);
% A(1:3,4:6)=eye(3)*dt;
% A(1:3,7:9)=eye(3)*dt^2/2;
% A(4:6,7:9)=eye(3)*dt;
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

function [A,Q]=transModelAcc(x,dt,palpha,pbeta,Da)
%数据准备
n=length(x);
A=eye(n,n);

%基本参数
v=x(4:6);   %速度分量
v2=v'*v;
if(v2<0.0001)
    v2=0.0001;
end
va=sqrt(v2);

exp_at=exp(-palpha*dt);
exp_at1=(1-exp_at)/palpha;
exp_at2=(dt-exp_at1)/palpha;

A(1:3,4:6)=eye(3)*(exp_at1-exp_at2*pbeta/va);
A(1:3,7:9)=eye(3)*exp_at2/va;
A(4:6,4:6)=eye(3)*(exp_at-exp_at1*pbeta/va);
A(4:6,7:9)=eye(3)*exp_at1/va;


% A=eye(n,n);
% A(1:3,4:6)=eye(3)*dt;
% A(1:3,7:9)=eye(3)*dt^2/2;
% A(4:6,7:9)=eye(3)*dt;

% 方差
Q=zeros(n,n);
Q(1:3,1:3)=Da*dt*exp_at2^2/5;
Q(1:3,4:6)=Da*dt*exp_at2*exp_at1/4;
Q(1:3,7:9)=Da*dt*exp_at2/5;

Q(4:6,4:6)=Da*dt*exp_at1^2/3;
Q(4:6,7:9)=Da*dt*exp_at1/2;
Q(7:9,7:9)=Da*dt;

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

% 5.计算新的轨迹
function Xs_up=MAPestimationPP(As,Qs,H,cut_t,Xs,Rzs_kp,Rzs_br,Zxp,Zbr,nxt)
%  trajectory MAP from observations
%  输入：离散状态转移矩阵集合As,状态转移协方差矩阵集合Qs,
%       观测协方差矩阵Rzs,观测矩阵H,截断时刻cut_t,上一时刻轨迹估计Xs,观测Zs
%  输出：状态轨迹X_re

n=length(As);  %状态时刻数
if(n==0)
    disp('error->in maxLikelihoodFilter 状态转移矩阵数目小于1 无法求解');
    return;
end

%只有少许状态，直接计算初始值
% if n<=5
%     Xs_up=0.1*(rand(10,n)-0.5);
%     for i=1:n
%         Xs_up(1:3,i)=positionFromObserv(Zs(:,i,:),Rzs(i,:));
%         Xs_up(7,i)=abs(Xs_up(7,i))*10;
%         if i>=2
%             Xs_up(4:6,i-1)=Xs_up(4:6,i-1)+(Xs_up(1:3,i)-Xs_up(1:3,i-1))/(time(i)-time(i-1));%速度
%         end
%     end
%     
%     if n>1
%         Xs_up(4:6,n)=Xs_up(4:6,n)+Xs_up(4:6,n-1);
%     end
%     return;
% end

%1.方差计算
%(1)计算状态方差的逆
n_use=n-cut_t;

Qs_inv=cell(n_use,1);

for i=1:n_use
    t=cut_t+i;  %时间戳
    
    Qs_inv{i}=inv(Qs{t});   %状态转移协方差的逆矩阵
end

if(cut_t>0) %前一个协方差阵求逆
    Qp_inv=inv(Qs{cut_t});
end

%(2)计算观测方差的逆
[nkp,~]=size(Rzs_kp);%关键点数
Rkp_inv=cell(nkp,1);
for i=1:nkp
    Rkp_inv{i}=inv(Rzs_kp{i});   %关键点逆矩阵
end

nbr=length(Rzs_br);
Rbr_inv=cell(nbr,1);
for i=1:nbr
    Rbr_inv{i}=inv(Rzs_br{i});   %障碍物逆矩阵
end

%2.系数矩阵与向量的各个分量
% (1)系数矩阵与向量的各个分量
m=size(Xs,1);

Mc=cell(n_use,1);   %M的对角元
Ms=cell(n_use,1);   %M的非对角元
b=zeros(m,n_use);   %每列同一个时刻，行方向为时间轴

% 第一个时刻
t=cut_t+1;
[Mc_b,bb]=barrierPara(Rbr_inv,Zbr,H,Xs(:,t));

xm_Mat=As{t}'*Qs_inv{1};
Ms{1}=xm_Mat;

b_cr=bb;
cMc=xm_Mat*As{t}+Mc_b;

if(cut_t>0)    %截断时刻生效
    cMc=cMc+Qp_inv; %t==cut_t
    b_cr=b_cr+Qp_inv*As{cut_t}*Xs(:,cut_t);
end

Mc{1}=cMc;
b(:,1)=b_cr;

% 后续时刻
for i=2:n_use-1
    t=cut_t+i;  %时间戳
    [Mc_b,bb]=barrierPara(Rbr_inv,Zbr,H,Xs(:,t));
    
    xm_Mat=As{t}'*Qs_inv{i};
    Ms{i}=xm_Mat;
    
    Mc{i}=xm_Mat*As{t}+Qs_inv{i-1}+Mc_b; %n_use<t<k
    b(:,i)=bb;
end

%t==k
t=cut_t+n_use;
[Mc_b,bb]=barrierPara(Rbr_inv,Zbr,H,Xs(:,t));

Mc{n_use}=Qs_inv{n_use-1}+Mc_b; %n_use=k
b(:,n_use)=bb;

%添加观测量
for j=1:nkp
    t=nxt(j);
    i=t-cut_t;
    
    mM=H'*Rkp_inv{j};
    Mc{i}=Mc{i}+mM*H;
    b(:,i)=b(:,i)+mM*Zxp(:,j);
end

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
    X_up(4:6,n_use)=(X_up(1:3,n_use)-Xs(1:3,t-2));
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

function [Mc_b,bb]=barrierPara(invRbr,Zbr,H,x)
%计算障碍物引起的方差
nbr=length(invRbr);
nx=length(x);
Mc_b=zeros(nx,nx);
bb=zeros(nx,1);
p=x(1:3);

for i=1:nbr
    pb=Zbr(:,i);
    dp=p-pb;
    
    %判断合理范围
    nrm2=dp'*invRbr{i}*dp;  %误差平方
    if nrm2>1.2
        continue;
    end
    
    %计算系数
    ex_nr=exp(-nrm2^4/5);
    ki=ex_nr/(1-ex_nr/1.02)*nrm2^3*H'*invRbr{i};
    
    Mc_b=Mc_b-ki*H;
    bb=bb-ki*pb;
end

end

% 计算状态转移误差
function dX_re=stransError(As,Xs,dX,cut_t)
n=length(As);
% 状态转移误差
dX_re=dX;
for i=cut_t+1:n-1
    dX_re(:,i)=As{i}*Xs(:,i)-Xs(:,i+1);
end
end

% 状态转移方差更新函数
function Qs=transCovUpdate(Qs,delta_x,cut_t,dt)
% 输入：状态转移矩阵Qs 当前轨迹状态转移误差delta_x 采样时刻time 
%       截断时刻cut_t 当前机动tm

nt=length(Qs);
if(nt<2)
    return;
end

time_use=nt-(cut_t+1);

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
ang=0.01;        %误差限

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

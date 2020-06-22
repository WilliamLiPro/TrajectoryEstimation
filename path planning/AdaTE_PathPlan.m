%���� adaptive trajectory estimation �Ĺ켣���ƺ���
%ʹ��DPTģ��
%����·���滮
function [As_up,Qs_up,cut_t,Xs_up,dX_up,preX]=AdaTE_PathPlan(As,Qs,Rzs_kp,Rzs_br,H,cut_t,Xs,dX,Zkp,Zbr,sitar,time,dt)
%@As ״̬ת�ƾ���
%@Qs �۲�Э����
%@Rzs_kp �ؼ��������ƫ�뷽��
%@Rzs_br �ϰ���ĳߴ��Ӧ����
%@H �۲����
%@cut_t �ض�ʱ��(��ǰʱ��)
%@Xs �켣״̬
%@dX ת�����
%@Zkp ·���ؼ�������
%@Zbr ·���ϰ�������
%@sitar ����
%@time ʱ��(�������key point ʱ��)
%@dt �켣����ʱ����

% ����׼��
palpha=sitar.alpha;
pbeta=sitar.beta;
Da=sitar.Da;
Dt=sitar.Dt;

[~,nx]=size(Xs);
[~,nzk]=size(Zkp);  %key point����

if(nzk<=1)
    preX=-1;    %�۲���Ϣ�������˲���������
    disp('error:�۲���Ϣ���㣬�޷�����˲�������۲���Ϣ');
    return;
end

%ѡȡ����2���ؼ���
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

%״̬��ʼ��
lx=size(Xs(:,nxt(1)),1);
if sum(Xs(:,nxt(1))==0)    %��ʼ������
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

if sum(Xs(:,nxt(2))==0)    %��ʼ������
    cur_x=Xs(:,nxt(1));
    cur_x(1:3)=Zxp(:,2);
    Xs(:,nxt(2))=cur_x+[0;0;0;rand(lx-3,1)*0.1-0.05];
end

% 1.����״̬ת�ƾ���
As_up=[As(1:nx);cell(nxt(2)-nx)];
Qs_up=[Qs(1:nx);cell(nxt(2)-nx)];
for i=cut_t+1:nxt(2)-1
    As_up{i}=AccModel(Xs(:,i),dt,palpha,pbeta);
%     As_up{i}=MotionModal2(Xs(:,i),dt,palpha,pbeta);
    preX=As_up{i}*Xs(:,i);  %Ԥ��
    
    if sum(Xs(:,i+1)==0)
        Xs(:,i+1)=preX+[0;0;0;rand(lx-3,1)*0.01-0.005];
    end
end

% 2.�������һ��Э����
if(cut_t==0||cut_t==time(nxp(1)-1))    %�յ����µĵ�
    for i=cut_t+1:nxt(2)-1
%         [As_up{i},Qs_up{i}]=transModel(Xs(:,i),dt,palpha,pbeta,Da,Dt);
        [As_up{i},Qs_up{i}]=transModelAcc(Xs(:,i),dt,palpha,pbeta,Da);
        
        preX=As_up{i}*Xs(:,i);  %Ԥ��
        
        if sum(Xs(:,i+1)==0)
            Xs(:,i+1)=preX+[0;0;0;rand(lx-3,1)*0.01-0.005];
        end
    end
end

% 4.���¹۲ⷽ��
%Rzs_up=observCovUpdate(Rzs_ob,Rzs,H,cut_t,[Xs,preX],Zs);

% 5.�����µĹ켣

if cut_t==0
    Xs_up=Xs;%��ʼ����Լ�������޷���⣬ֱ���˳�
else
    Xs_up=MAPestimationPP(As_up,Qs_up,H,cut_t,Xs,Rzs_kp(nxp(1):nxp(2)),Rzs_br,Zxp,Zbr,nxt);
end

% 6.����״̬ת�Ʒ���
dX_up=stransError(As_up,Xs_up,dX,cut_t);  %����ת�����
Qs_up=transCovUpdate(Qs_up,dX_up,cut_t,dt);

% 7.�ض�ʱ�̸���
% Fiy=diag([1,1,1,1/4,1/4,1/4,1/9,1/9,1/9,1/9]);
% cut_t=cutTimeUpdate([Xs,preX],Xs_up,Fiy,cut_t);
end

function A=AccModel(x,dt,palpha,pbeta)
%  ����״̬ת�ƾ���A
%  ���룺 ��һʱ��״̬����x ʱ����dt �������� palpha,pbeta
%  �����״̬ת�ƾ���A

%����׼��
n=length(x);
A=eye(n,n);

%��������
v=x(4:6);   %�ٶȷ���
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
%  ����״̬ת�ƾ���A,ʹ��ƫ�������Ƽ���
%  ���룺 ״̬����x,xp ��һ���ٶ�vn ʱ����dt �������� palpha,pbeta
%  �����״̬ת�ƾ���A

%����׼��
n=length(x);
A=eye(n,n);

%�ٶ�
v=x(4:6);   %�ٶȷ���
va=sqrt(v'*v);
if va<0.00001
    va=0.00001;
end

%����
w=x(7);
if w<-pbeta^2/(4*palpha);
    w=0.000001-pbeta^2/(4*palpha);
end

%������һʱ���ٶ�
sqt_baw=sqrt(pbeta^2+4*palpha*w);
r=pbeta/sqt_baw;
kesy1=(pbeta+sqt_baw)/(2*palpha);
kesy2=-2*w/(pbeta+sqt_baw);

va_n=va;
f2=(1+r)*log(abs(va+kesy1))+(1-r)*log(abs(va+kesy2))-2*palpha*dt;  %ǰһʱ��ϵ��
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

%��������
exp_ra=exp(-2*palpha*dt/(1-r));
exp_ra1=(1-exp_ra)/(2*palpha/(1-r));

hk=abs((va_n+kesy1)/(va+kesy1))^(-(1+r)/(1-r));
gk=2/(pbeta+sqt_baw)*(1-exp_ra*hk);
gk1=2/(pbeta+sqt_baw)*(dt-exp_ra1*hk);

%�������
vx=zeros(3,3);%�������
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
%����׼��
n=length(x);
A=eye(n,n);

%�ٶ�
v=x(4:6);   %�ٶȷ���
va=sqrt(v'*v);
if va<0.00001
    va=0.00001;
end

%����
w=x(7);
if w<-pbeta^2/(4*palpha);
    w=0.000001-pbeta^2/(4*palpha);
end

%������һʱ���ٶ�
sqt_baw=sqrt(pbeta^2+4*palpha*w);
r=pbeta/sqt_baw;
kesy1=(pbeta+sqt_baw)/(2*palpha);
kesy2=-2*w/(pbeta+sqt_baw);

va_n=va;
f2=(1+r)*log(abs(va+kesy1))+(1-r)*log(abs(va+kesy2))-2*palpha*dt;  %ǰһʱ��ϵ��
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

%��������
exp_ra=exp(-2*palpha*dt/(1-r));
exp_ra1=(1-exp_ra)/(2*palpha/(1-r));

hk=abs((va_n+kesy1)/(va+kesy1))^(-(1+r)/(1-r));
gk=2/(pbeta+sqt_baw)*(1-exp_ra*hk);
gk1=2/(pbeta+sqt_baw)*(dt-exp_ra1*hk);

%�������
vx=zeros(3,3);%�������
vx(1,2)=-v(3); vx(1,3)=v(2);
vx(2,1)=v(3);  vx(2,3)=-v(1);
vx(3,1)=-v(2); vx(3,2)=v(1);

A(1:3,4:6)=eye(3)*(exp_ra1*hk);
A(1:3,7)=gk1*v/va;
A(1:3,8:10)=-(exp_ra1*hk+gk1*w/va)*dt*vx;
A(4:6,4:6)=eye(3)*(exp_ra*hk);
A(4:6,7)=gk*v/va;
A(4:6,8:10)=-(exp_ra*hk+gk*w/va)*dt*vx;

% 4.Э����
Q=zeros(n,n);

%��˾���
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

%�Գƻ�
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
%����׼��
n=length(x);
A=eye(n,n);

%��������
v=x(4:6);   %�ٶȷ���
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

% ����
Q=zeros(n,n);
Q(1:3,1:3)=Da*dt*exp_at2^2/5;
Q(1:3,4:6)=Da*dt*exp_at2*exp_at1/4;
Q(1:3,7:9)=Da*dt*exp_at2/5;

Q(4:6,4:6)=Da*dt*exp_at1^2/3;
Q(4:6,7:9)=Da*dt*exp_at1/2;
Q(7:9,7:9)=Da*dt;

%�Գƻ�
for i=1:n
    for j=i+1:n
        Q(j,i)=Q(i,j);
    end
end

for i=1:n
    Q(i,i)=Q(i,i)*n/2;
end

end

% 5.�����µĹ켣
function Xs_up=MAPestimationPP(As,Qs,H,cut_t,Xs,Rzs_kp,Rzs_br,Zxp,Zbr,nxt)
%  trajectory MAP from observations
%  ���룺��ɢ״̬ת�ƾ��󼯺�As,״̬ת��Э������󼯺�Qs,
%       �۲�Э�������Rzs,�۲����H,�ض�ʱ��cut_t,��һʱ�̹켣����Xs,�۲�Zs
%  �����״̬�켣X_re

n=length(As);  %״̬ʱ����
if(n==0)
    disp('error->in maxLikelihoodFilter ״̬ת�ƾ�����ĿС��1 �޷����');
    return;
end

%ֻ������״̬��ֱ�Ӽ����ʼֵ
% if n<=5
%     Xs_up=0.1*(rand(10,n)-0.5);
%     for i=1:n
%         Xs_up(1:3,i)=positionFromObserv(Zs(:,i,:),Rzs(i,:));
%         Xs_up(7,i)=abs(Xs_up(7,i))*10;
%         if i>=2
%             Xs_up(4:6,i-1)=Xs_up(4:6,i-1)+(Xs_up(1:3,i)-Xs_up(1:3,i-1))/(time(i)-time(i-1));%�ٶ�
%         end
%     end
%     
%     if n>1
%         Xs_up(4:6,n)=Xs_up(4:6,n)+Xs_up(4:6,n-1);
%     end
%     return;
% end

%1.�������
%(1)����״̬�������
n_use=n-cut_t;

Qs_inv=cell(n_use,1);

for i=1:n_use-1
    t=cut_t+i;  %ʱ���
    
    Qs_inv{i}=inv(Qs{t})+eye(9)*0.00000001;   %״̬ת��Э����������
end

if(cut_t>0) %ǰһ��Э����������
    Qp_inv=inv(Qs{cut_t})+eye(9)*0.00000001;
end

%(2)����۲ⷽ�����
[nkp,~]=size(Rzs_kp);%�ؼ�����
Rkp_inv=cell(nkp,1);
for i=1:nkp
    Rkp_inv{i}=inv(Rzs_kp{i});   %�ؼ��������
end

nbr=length(Rzs_br);
Rbr_inv=cell(nbr,1);
for i=1:nbr
    Rbr_inv{i}=inv(Rzs_br{i});   %�ϰ��������
end

%2.ϵ�������������ĸ�������
% (1)ϵ�������������ĸ�������
m=size(Xs,1);

Mc=cell(n_use,1);   %M�ĶԽ�Ԫ
Ms=cell(n_use,1);   %M�ķǶԽ�Ԫ
b=zeros(m,n_use);   %ÿ��ͬһ��ʱ�̣��з���Ϊʱ����

% ��һ��ʱ��
t=cut_t+1;
[Mc_b,bb]=barrierPara(Rbr_inv,Zbr,H,Xs(:,t));

xm_Mat=As{t}'*Qs_inv{1};
Ms{1}=xm_Mat;

b_cr=bb;
cMc=xm_Mat*As{t}+Mc_b;

if(cut_t>0)    %�ض�ʱ����Ч
    cMc=cMc+Qp_inv; %t==cut_t
    b_cr=b_cr+Qp_inv*As{cut_t}*Xs(:,cut_t);
end

Mc{1}=cMc;
b(:,1)=b_cr;

% ����ʱ��
for i=2:n_use-1
    t=cut_t+i;  %ʱ���
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

%��ӹ۲���
for j=1:nkp
    t=nxt(j);
    i=t-cut_t;
    
    mM=H'*Rkp_inv{j};
    Mc{i}=Mc{i}+mM*H;
    b(:,i)=b(:,i)+mM*Zxp(:,j);
end

%2.���ݹ۲�ֵУ���켣
% (1)������Ԫ
for i=1:n_use-1
    cg=Ms{i}'/Mc{i};
    Mc{i+1}=Mc{i+1}-cg*Ms{i};
    b(:,i+1)=b(:,i+1)+cg*b(:,i);
end

% (3)��������
X_up=zeros(m,n_use);	%X���²���
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

% �������
Xs_up=[Xs(:,1:cut_t),X_up];
end

function [Mc_b,bb]=barrierPara(invRbr,Zbr,H,x)
%�����ϰ�������ķ���
nbr=length(invRbr);
nx=length(x);
Mc_b=zeros(nx,nx);
bb=zeros(nx,1);
p=x(1:3);

for i=1:nbr
    pb=Zbr(:,i);
    dp=p-pb;
    
    %�жϺ���Χ
    nrm2=dp'*invRbr{i}*dp;  %���ƽ��
    if nrm2>1.5
        continue;
    end
    
    %����ϵ��
    ex_nr=exp(-nrm2^4/5);
    ki=ex_nr/(1-ex_nr/1.02)*nrm2^3*H'*invRbr{i};
    
    Mc_b=Mc_b-ki*H;
    bb=bb-ki*pb;
end

end

% ����״̬ת�����
function dX_re=stransError(As,Xs,dX,cut_t)
n=length(As);
% ״̬ת�����
dX_re=dX;
for i=cut_t+1:n-1
    dX_re(:,i)=As{i}*Xs(:,i)-Xs(:,i+1);
end
end

% ״̬ת�Ʒ�����º���
function Qs=transCovUpdate(Qs,delta_x,cut_t,dt)
% ���룺״̬ת�ƾ���Qs ��ǰ�켣״̬ת�����delta_x ����ʱ��time 
%       �ض�ʱ��cut_t ��ǰ����tm

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

% �ض�ʱ�̸���
function cut_t=cutTimeUpdate(Xs_pre,Xs_ev,Fiy,cut_tp)
%  ���룺״̬�켣Ԥ��Xs_pre ״̬�켣����Xs_ev Ȩ�ؾ���Fiy ǰһ���ض�ʱ��cut_tp
%  ������ض�ʱ��cut_t

[~,n]=size(Xs_pre);
ang=0.01;        %�����

sum_eer=0;

for i=cut_tp+1:2:n-40   %������ʱ��һ�β���
    error=Xs_pre(:,i)-Xs_ev(:,i);
    if sum_eer==0
        sum_eer=error'*Fiy*error;
    else
        sum_eer=0.95*sum_eer+0.05*error'*Fiy*error;
    end
    if(sum_eer>ang)    %ƽ����������ֵ
        cut_t=max(i-20,0);
        return;
    end
end
cut_t=max(n-50,0);
end

function p=positionFromObserv(zs,Rz)
lz=length(Rz);%�۲�ʱ������ͬһʱ�̹۲����

p=zeros(3,1);
sum_w=zeros(3,3);

for j=1:lz
    iv_Rz=inv(Rz{j});
    p=p+iv_Rz*zs(:,j);
    sum_w=sum_w+iv_Rz;
end

p=sum_w\p;
end

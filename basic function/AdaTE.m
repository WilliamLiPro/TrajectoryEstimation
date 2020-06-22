%���� adaptive trajectory estimation �Ĺ켣���ƺ���
%ʹ��DPTģ��
function [As_up,Qs_up,Rzs_up,cut_t,Xs_up,dX_up,preX]=AdaTE(As,Qs,Rzs_ob,Rzs,H,cut_t,Xs,dX,Zs,sitar,time)
%@As ״̬ת�ƾ���
%@Qs �۲�Э����
%@Rzs_ob �۲ⷽ��ĳ�ʼֵ
%@Rzs ������Ĺ۲ⷽ��
%@H �۲����
%@cut_t �ض�ʱ��
%@Xs �켣״̬
%@dX ת�����
%@Zs �۲�
%@sitar ����
%@time ʱ��

% ����׼��
palpha=sitar.alpha;
pbeta=sitar.beta;
Da=sitar.Da;
Dt=sitar.Dt;

[~,nx]=size(Xs);
[~,nz]=size(Zs);

if(nz<=nx)
    preX=-1;    %�۲���Ϣ�������˲���������
    disp('error:�۲���Ϣ���㣬�޷�����˲�������۲���Ϣ');
    return;
end

% 1.����״̬ת�ƾ���
As_up=[As(1:nx);cell(1)];
Qs_up=[Qs(1:nx);cell(1)];
for i=cut_t+1:nx-1
    dt=time(i+1)-time(i);
    As_up{i}=MotionModal2(Xs(:,i),dt,palpha,pbeta);
end

% 2.�������һ��Э����
if(nx>0)
    dt=time(nx+1)-time(nx);
    [As_up{nx},Qs_up{nx}]=transModel(Xs(:,nx),dt,palpha,pbeta,Da,Dt);
end

% 3.�켣Ԥ��
if(nx==0)
    preX=zeros(10,1);
    preX(1:3,1)=Zs(1:3,1);
    preX(4:6,1)=(rand(3,1)-0.5);
else
    preX=As_up{nx}*Xs(:,nx);
end

% 4.���¹۲ⷽ��
Rzs_up=observCovUpdate(Rzs_ob,Rzs,H,cut_t,[Xs,preX],Zs);

% 5.�����µĹ켣
Xs_up=MAPestimation(As_up,Qs_up,Rzs_up,H,cut_t,Xs,Zs,time);

% 6.����״̬ת�Ʒ���
dX_up=stransError(As_up,Xs_up,dX,cut_t);  %����ת�����
Qs_up=transCovUpdate(Qs_up,dX_up,time,cut_t);

% 7.�ض�ʱ�̸���
Fiy=diag([1,1,1,1/4,1/4,1/4,1/9,1/9,1/9,1/9]);
cut_t=cutTimeUpdate([Xs,preX],Xs_up,Fiy,cut_t);

end

function A=MotionModal(x,dt,palpha,pbeta)
%  ����״̬ת�ƾ���A,ƫת����F
%  ���룺 ��һʱ��״̬����x ʱ����dt �������� palpha,pbeta
%  �����״̬ת�ƾ���A

n=length(x);
A=eye(n,n);

w=x(7);
if w<0
    w=0;
end

%��������
v=x(4:6);   %�ٶȷ���
v2=v'*v;
if(v2<0.000001)
    v2=0.000001;
end
va=sqrt(v2);

% 1.����״̬ת����ƫ�ƾ���
if w<pbeta^2/(palpha*4)
    r=pbeta/sqrt(pbeta^2+4*palpha*w);
    r=2*palpha/(1+r);
    
    exp_r=exp(-r*dt);%ָ����
    exp_r1=(1-exp_r)/r;%ָ����һ�׻���
    texp_r=dt*exp_r;%ʱ��ָ����
    texp_r1=(exp_r1-texp_r)/r;%ʱ��ָ����
    
    peta_v=pbeta/va;
    
    %λ��
    for i=1:3
        A(i,i+3)=exp_r1-dt^2*peta_v/2;
    end
    A(1:3,8:10)=dt^3*peta_v/3-texp_r1;
    
    % �ٶ�����
    for i=4:6
        A(i,i)=exp_r-dt*peta_v;
    end
    A(4:6,8:10)=dt^2*peta_v-texp_r;
else
    exp_a=exp(-palpha*dt);    %ָ����
    exp_a1=(1-exp_a)/palpha;%ָ����һ�׻���
    texp_a=dt*exp_a;    %ʱ��ָ����
    texp_a1=(exp_a1-texp_a)/palpha;%ʱ��ָ��һ�׻���
    
    para_gt=v2+dt*w/2;
    para_w=v/para_gt;       %���ʷ���
    para_ct=zeros(3,3);%�������
    para_ct(1,2)=-v(3); para_ct(1,3)=v(2);
    para_ct(2,1)=v(3);  para_ct(2,3)=-v(1);
    para_ct(3,1)=-v(2); para_ct(3,2)=v(1);
    
    para_ct=para_ct*(1+dt*w/para_gt);
    
    for i=1:3
        A(i,i+3)=exp_a1;
    end
    
    A(1:3,7)=texp_a1*para_w;
    A(1:3,8:10)=-texp_a1*para_ct;
    
    % �ٶ�����
    for i=4:6
        A(i,i)=exp_a;
    end
    
    A(4:6,7)=texp_a*para_w;
    A(4:6,8:10)=-texp_a*para_ct;
end
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

function A=MotionModal3(x,dt,palpha,pbeta)
%  ����״̬ת�ƾ���A,ʹ��ƫ�������Ƽ���
%  ���룺 ״̬����x,xp ��һ���ٶ�vn ʱ����dt �������� palpha,pbeta
%  �����״̬ת�ƾ���A

%1.����׼��
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
if abs(w)<0.01
    w=0.01*w/abs(w);
end

%2.������һʱ���ٶ�
ra=w/pbeta;

va_n=va;
f2=va+ra*log(abs(va-ra))-pbeta*dt;  %ǰһʱ��ϵ��
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

%3.״̬ת�ƾ���
%��������
exp_rb=exp(-pbeta*dt/ra);%ָ����
exp_rb1=(1-exp_rb)*ra/pbeta;%ָ����һ�׻���

hk=exp((va-va_n)/ra);
gk=(1-exp_rb*hk)/(pbeta*va);
gk1=(dt-exp_rb1*hk)/(pbeta*va);

%����
vx=zeros(3,3);%�������
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

function [A,Q]=transModel2(x,dt,palpha,pbeta,Da,Dt)
%1.����׼��
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

%2.������һʱ���ٶ�
ra=w/pbeta;

va_n=va;
f2=va+ra*log(abs(va-ra))-pbeta*dt;  %ǰһʱ��ϵ��
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

%3.״̬ת�ƾ���
%��������
exp_rb=exp(-pbeta*dt/ra);%ָ����
exp_rb1=(1-exp_rb)*ra/pbeta;%ָ����һ�׻���

hk=exp((va-va_n)/ra);
gk=(1-exp_rb*hk)/(pbeta*va);
gk1=(dt-exp_rb1*hk)/(pbeta*va);

%����
vx=zeros(3,3);%�������
vx(1,2)=-v(3); vx(1,3)=v(2);
vx(2,1)=v(3);  vx(2,3)=-v(1);
vx(3,1)=-v(2); vx(3,2)=v(1);

A(1:3,4:6)=(exp_rb1*hk)*eye(3);
A(1:3,7)=gk1*v;
A(1:3,8:10)=-(exp_rb1*hk+gk1*w)*dt*vx;

A(4:6,4:6)=(exp_rb*hk)*eye(3);
A(4:6,7)=gk*v;
A(4:6,8:10)=-(exp_rb*hk+gk*w)*dt*vx;

% 4.Э����
Q=zeros(n,n);

%��˾���
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

% 4.���¹۲ⷽ��
function Rzs=observCovUpdate(Rzs_ob,Rzs,H,cut_t,Xs,Zs)
[nr,lz]=size(Rzs);
mz=size(Zs,1);

if(nr<1)
    return;
end

% ����۲����
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

% 5.�����µĹ켣
function Xs_up=MAPestimation(As,Qs,Rzs,H,cut_t,Xs,Zs,time)
%  trajectory MAP from observations
%  ���룺��ɢ״̬ת�ƾ��󼯺�As,״̬ת��Э������󼯺�Qs,
%       �۲�Э�������Rzs,�۲����H,�ض�ʱ��cut_t,��һʱ�̹켣����Xs,�۲�Zs
%  �����״̬�켣X_re

n=length(As);  %״̬ʱ����
[nz,lz]=size(Rzs);%�۲�ʱ������ͬһʱ�̹۲����
if(nz<n)
    disp('error->in maxLikelihoodFilter �۲�������٣��޷����');
    return;
end
if(n==0)
    disp('error->in maxLikelihoodFilter ״̬ת�ƾ�����ĿС��1 �޷����');
    return;
end

%ֻ������״̬��ֱ�Ӽ����ʼֵ
if n<=5
    Xs_up=0.1*(rand(10,n)-0.5);
    for i=1:n
        Xs_up(1:3,i)=positionFromObserv(Zs(:,i,:),Rzs(i,:));
        Xs_up(7,i)=abs(Xs_up(7,i))*10;
        if i>=2
            Xs_up(4:6,i-1)=Xs_up(4:6,i-1)+(Xs_up(1:3,i)-Xs_up(1:3,i-1))/(time(i)-time(i-1));%�ٶ�
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

%1.�������״̬�ĵ������
%(1)���㷽�����
n_use=n-cut_t;

Qs_inv=cell(n_use,1);
Rs_inv=cell(n_use,lz);

for i=1:n_use
    t=cut_t+i;  %ʱ���
    
    Q=Qs{t};    %״̬ת��Э����
    if i<n_use
        Qs_inv{i}=inv(Q)+eye(10)*0.00000001;   %�����
    end
    
    for j=1:lz
        Rs_inv{i,j}=inv(Rzs{t,j});  %�۲�Э������
    end
end

if(cut_t>0) %ǰһ��Э����������
    Qp_inv=inv(Qs{cut_t});
end

% (2)ϵ�������������ĸ�������
m=size(Xs,1);

Mc=cell(n_use,1);   %M�ĶԽ�Ԫ
Ms=cell(n_use,1);   %M�ķǶԽ�Ԫ
b=zeros(m,n_use);   %ÿ��ͬһ��ʱ�̣��з���Ϊʱ����

% ��һ��ʱ��
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

if(cut_t>0)    %�ض�ʱ����Ч
    Mc{1}=Mc{1}+Qp_inv; %t==cut_t
    b(:,1)=b(:,1)+Qp_inv*As{cut_t}*Xs(:,cut_t);
end

% ����ʱ��
for i=2:n_use-1
    t=cut_t+i;  %ʱ���
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

% �������
Xs_up=[Xs(:,1:cut_t),X_up];
end

% ����״̬ת�����
function dX_re=stransError(As,Xs,dX,cut_t)
n=length(As);
% ״̬ת�����
dX_re=[dX,zeros(10,1)];
for i=cut_t+1:n-1
    dX_re(:,i)=As{i}*Xs(:,i)-Xs(:,i+1);
end
end

% ״̬ת�Ʒ�����º���
function Qs=transCovUpdate(Qs,delta_x,time,cut_t)
% ���룺״̬ת�ƾ���Qs ��ǰ�켣״̬ת�����delta_x ����ʱ��time 
%       �ض�ʱ��cut_t ��ǰ����tm

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

% �ض�ʱ�̸���
function cut_t=cutTimeUpdate(Xs_pre,Xs_ev,Fiy,cut_tp)
%  ���룺״̬�켣Ԥ��Xs_pre ״̬�켣����Xs_ev Ȩ�ؾ���Fiy ǰһ���ض�ʱ��cut_tp
%  ������ض�ʱ��cut_t

[~,n]=size(Xs_pre);
ang=0.000001;        %�����

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

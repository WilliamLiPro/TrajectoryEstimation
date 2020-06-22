%���� adaptive trajectory estimation �Ĺ켣���ƺ���
%ʹ�ü��ٶ�ģ��
function [As_up,Qs_up,Rzs_up,cut_t,Xs_up,dX_up,preX]=AdaTEAcc(As,Qs,Rzs_ob,Rzs,H,cut_t,Xs,dX,Zs,sitar,time)
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
    As_up{i}=AccModel(Xs(:,i),dt,palpha,pbeta);
end

% 2.�������һ��Э����
if(nx>0)
    dt=time(nx+1)-time(nx);
    [As_up{nx},Qs_up{nx}]=transModel(Xs(:,nx),dt,palpha,pbeta,Da);
end

% 3.�켣Ԥ��
if(nx==0)
    preX=zeros(9,1);
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
% Qs_up=transCovUpdate(Qs_up,dX_up,time,cut_t);

% 7.�ض�ʱ�̸���
Fiy=diag([1,1,1,1/4,1/4,1/4,1/9,1/9,1/9]);
cut_t=cutTimeUpdate([Xs,preX],Xs_up,Fiy,cut_t);

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

end

function [A,Q]=transModel(x,dt,palpha,pbeta,Da)
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
if n<=6
    Xs_up=0.1*(rand(9,n)-0.5);
    for i=1:n
        Xs_up(1:3,i)=positionFromObserv(Zs(:,i,:),Rzs(i,:));
        
        if i>=2
            Xs_up(4:6,i-1)=Xs_up(4:6,i-1)+(Xs_up(1:3,i)-Xs_up(1:3,i-1))/(time(i)-time(i-1));%�ٶ�
        end
    end
    
    if n>=2
        Xs_up(4:6,i)=Xs_up(4:6,i)+Xs_up(4:6,i-1);
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
    
    if i<n_use
        Q=Qs{t};    %״̬ת��Э����
        Qs_inv{i}=inv(Q)+eye(9)*0.00000001;   %�����
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
b_cr=Rs_inv{1,1}*Zs(:,1,1);
sum_Rs=Rs_inv{1,1};
t=cut_t+1;
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
        X_up(:,1:i)=Xs(:,cut_t+1:cut_t+i)+0.04*(rand(9,i)-0.5);
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
dX_re=[dX,zeros(9,1)];
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
ang=0.0001;        %�����

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

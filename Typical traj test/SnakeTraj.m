%% snake like trajectory
%   �����ֲ��۲ⶪʧ
%   ���룺ƽ�й۲����nz ����Rz Ư��ʱ��ts te Ư�ƹ۲����rate
%   ������켣��ֵreal �۲�켣obs �۲�ʱ��time

function [real,obs,Rzs,time]=SnakeTraj(nz,Rz,ts,te,rate)
n=600;
time=[1:n];
real=zeros(3,n);
Rzs=cell(n,nz);

% 1.���ɹ켣��ֵ
% ������
for i=1:200
    real(1,i)=0.2*i+40*sin(i/31.83-pi/4);
    real(2,i)=0.2*i+40*cos(i/31.83-pi/4);
end

% ���ٰ�Բ
sx=real(1,200)+22*sqrt(2);sy=real(2,200)-22*sqrt(2);
for i=201:300
    s=i-200;
    real(1,i)=sx+(44+s/10)*sin(s/31.83-pi/4);
    real(2,i)=sy+(44+s/10)*sin(s/31.83+pi/4);
    real(3,i)=10*cos(s/31.83)-10;
end

% c����
v0=(real(:,300)-real(:,299))/(time(300)-time(299));
v=v0;
for i=301:400
    s=i-300;
    real(:,i)=real(:,i-1)+(time(i)-time(i-1))*v;
    v(1)=v0(1)*(1-s/50);
    v(2)=(101-s)/100*v0(2);
    v(3)=(2500-(i-350)^2)/5000;
end

% ˫Բ
v0=v;
r=1/(pi/50+1/50)*sqrt(sum(v0.^2));
sx=real(1,400);sy=real(2,400)+r;

for i=401:600
    s=i-400;
    real(1,i)=sx+r*(1+s/200)*sin(s/15.9155)+r*(1+s/600)*s/50;
    real(2,i)=sy-r*(1+s/200)*cos(s/15.9155);
    real(3,i)=real(3,400)*(3+cos(s/10))/4;
end

% 2.�����켣
% ���������켣
obs=zeros(3,n,nz);
root_Rz=sqrtm(Rz);
for i=1:n
    for j=1:nz
        obs(:,i,j)=real(:,i)+root_Rz*randn(3,1);
        Rzs{i,j}=Rz;
    end
end

% �۲ⶪʧ
Rz_lack=diag([1000000,1000000,1000000]);
for i=250:300
    for j=1:nz
        Rzs{i,j}=Rz_lack;
    end
end

% Ư��
if rate<=0
    return;
end

if(nz*rate<1)
    dt=fix(1/(nz*rate));
    for i=ts:dt:te
        obs(:,1,i)=obs(:,1,i)+[4;4;0]*sqrt(i-ts);
    end
else
    sift_n=fix(nz*rate);
    for i=ts:te
        for j=1:sift_n
            obs(:,i,j)=obs(:,i,j)+[4;4;0]*sqrt(i-ts);
        end
    end
end
end
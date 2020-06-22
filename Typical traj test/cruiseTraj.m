%% cruise trajectory
%   输入：平行观测个数nz 方差Rz 漂移时段ts te 漂移观测比例rate
%   输出：轨迹真值real 观测轨迹obs 观测时间time
function [real,obs,time]=cruiseTraj(nz,Rz,ts,te,rate)

n=500;
time=[1:500];
real=zeros(3,n);

% 1.生成轨迹真值
% 直线
for i=1:50
    real(1,i)=0;
    real(2,i)=50-i;
end

% 3/4圆
for i=51:200
    sita=(i-50)/32;
    real(1,i)=32*cos(sita)-32;
    real(2,i)=-32*sin(sita);
end

% 直线
for i=201:300
    v=1+(i-200)/100;
    real(1,i)=real(1,i-1)+v;
    real(2,i)=real(2,i-1);
end

% 3/4圆
sx=real(1,300);sy=real(2,300);
for i=301:450
    sita=(i-300)/32;
    real(1,i)=sx+sin(sita)*64;
    real(2,i)=sy+50-cos(sita)*50;
end

% 直线
sx=real(1,450);
for i=451:n
    v=50/32-(i-450)/100;
    real(1,i)=sx;
    real(2,i)=real(2,i-1)-v;
end

% 2.噪声轨迹
% 基础噪声轨迹
obs=zeros(3,n,nz);
root_Rz=sqrtm(Rz);
for i=1:n
    for j=1:nz
        obs(:,i,j)=real(:,i)+root_Rz*randn(3,1);
    end
end

% 漂移
if rate<=0
    return;
end

if(nz*rate<1)
    dt=fix(1/(nz*rate));
    for i=ts:dt:te
        obs(:,i,1)=obs(:,i,1)+[8;8;0]*sqrt(i-ts);
    end
else
    sift_n=fix(nz*rate);
    for i=ts:te
        for j=1:sift_n
            obs(:,i,j)=obs(:,i,j)+[8;8;0]*sqrt(i-ts);
        end
    end
end

end
%% S�͹켣���ɺ���
%   ���룺ƽ�й۲����nz ����Rz Ư��ʱ��ts te Ư�ƹ۲����rate
%   ������켣��ֵreal �۲�켣obs �۲�ʱ��time

function [real,obs,time]=swayingCurve(nz,Rz,ts,te,rate)
n=500;
time=[1:500];
real=zeros(3,n);

% 1.���ɹ켣��ֵ
% �켣ΪS�ͣ�����𽥼���,�ڶ���������
for i=1:n
    real(1,i)=i+i^2/n/8;
    real(2,i)=i/10*sin(i/20);
end

% 2.�����켣
% ���������켣
obs=zeros(3,n,nz);
root_Rz=sqrtm(Rz);
for i=1:n
    for j=1:nz
        obs(:,i,j)=real(:,i)+root_Rz*randn(3,1);
    end
end

% Ư��
if rate<=0
    return;
end

if(nz*rate<1)
    dt=fix(1/(nz*rate));
    for i=ts:dt:te
        obs(:,1,i)=obs(:,1,ts-1)+root_Rz*randn(3,1)+sqrt([i-ts;(i-ts);0])*8;
    end
else
    sift_n=fix(nz*rate);
    for i=ts:te
        for j=1:sift_n
            obs(:,i,j)=obs(:,ts-1,j)+root_Rz*randn(3,1)+sqrt([i-ts;(i-ts);0])*8;
        end
    end
end

% if(nz*rate<1)
%     dt=fix(1/(nz*rate));
%     for i=ts:dt:te
%         obs(:,i,1)=obs(:,i,1)+[0;5*sqrt(i-ts);0];
%     end
% else
%     sift_n=fix(nz*rate);
%     for i=ts:te
%         for j=1:sift_n
%             obs(:,i,j)=obs(:,i,j)+[0;10*sqrt(i-ts);0];
%         end
%     end
% end
end
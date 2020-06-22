%%  trjectory estimation experiment
%   OAO test
%   contrast algorithm: EKF-RTS-Acc,EKF-RTS-DPT,UKF-URTS-Acc,UKF-URTS-DPT

%   result£º
%   observation trajectory RMSE-time curve
%   table of error (mean RMSE in each group)
%   graph of error (distribution of mean RMSE)

%%  path
addpath('../basicFunction');

%%  cruise trajectory without colored noise
nz=4;
Rz=diag([256,259,64]);
re_cw=trajectoryTest('cruise',nz,Rz,200,150,0);

%%  cruise trajectory with colored noise
re_cc=trajectoryTest('cruise',nz,Rz,300,450,0.25);

%%  swaying curve without colored noise
nz=4;
Rz=diag([256,256,64]);
re_sw=trajectoryTest('swaying',nz,Rz,200,150,0);

%%  swaying curve with colored noise
re_sc=trajectoryTest('swaying',nz,Rz,200,450,0.25);

%%  snake-like trajectory without colored noise
nz=4;
Rz=diag([64,69,63]);
re_slw=trajectoryTest('snake',nz,Rz,200,150,0);

%%  snake-like trajectory with colored noise
re_slc=trajectoryTest('snake',nz,Rz,450,600,0.25);

%%  draw the bars of statistics
% the min-max of RMSE
e_min=zeros(6,6);
e_max=zeros(6,6);

e_min(:,1)=re_cw.error(:,2);
e_min(:,2)=re_cc.error(:,2);
e_min(:,3)=re_sw.error(:,2);
e_min(:,4)=re_sc.error(:,2);
e_min(:,5)=re_slw.error(:,2);
e_min(:,6)=re_slc.error(:,2);

e_max(:,1)=re_cw.error(:,3);
e_max(:,2)=re_cc.error(:,3);
e_max(:,3)=re_sw.error(:,3);
e_max(:,4)=re_sc.error(:,3);
e_max(:,5)=re_slw.error(:,3);
e_max(:,6)=re_slc.error(:,3);

e_std=zeros(6,6);
for i=1:6
    e_std(i,1)=std(re_cw.RMSE{i});
    e_std(i,2)=std(re_cc.RMSE{i});
    e_std(i,3)=std(re_sw.RMSE{i});
    e_std(i,4)=std(re_sc.RMSE{i});
    e_std(i,5)=std(re_slw.RMSE{i});
    e_std(i,6)=std(re_slc.RMSE{i});
end


figure(10);
x=zeros(6,6);
for i=1:6
    for j=1:6
        x(i,j)=i+0.13*(j-6);
    end
end
y=zeros(6,6);
y(:,1)=re_cw.error(:,1);
y(:,2)=re_cc.error(:,1);
y(:,3)=re_sw.error(:,1);
y(:,4)=re_sc.error(:,1);
y(:,5)=re_slw.error(:,1);
y(:,6)=re_slc.error(:,1);
y(y>0.4)=0.4;

bar(x,y');hold on;
errorbar(x,y',e_std');

set(gca,'XTickLabel',{'a','b','c','d','e','f'});
ylabel('normalized RMSE');
axis([0.2,9,0,0.4]);
legend('AdaTE-PLS','MAP-CT','EKF-RTS-Acc','EKF-RTS-PLS','UKF-URTS-Acc','UKF-URTS-PLS');

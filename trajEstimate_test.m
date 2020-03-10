% test: trajectory estimation under noise observation
% Online adaptive trajectory estimation
% complexity < O(k): where k is the time

addpath('basic function','Typical traj test');

%% params
% trajectory type: trajType
% observation number at same time: nz
% observation covariance: Rz
% time range of colored noise: crn_st,crn_ft
% rate of colored noise: ratio
% number of observer: 4

trajType='snake';
nz=4;
Rz=diag([256,259,64]);
crn_st=450;
crn_ft=600;
ratio=0.25;

%% 1.create ground truth and observation
switch trajType
    case 'cruise'
        disp(['standard trajectory estimation 1: cruise trajectory']);
        [real,obs,time]=cruiseTraj(nz,Rz,crn_st,crn_ft,ratio);
    case 'swaying'
        disp(['standard trajectory estimation 2: swaying curve']);
        [real,obs,time]=swayingCurve(nz,Rz,crn_st,crn_ft,ratio);
    case 'snake'
        disp(['standard trajectory estimation 3: snake like trajectory']);
        [real,obs,rzs,time]=SnakeTraj(nz,Rz,crn_st,crn_ft,ratio);
end

x_min=min(real(1,:));x_max=max(real(1,:));
y_min=min(real(2,:));y_max=max(real(2,:));

x_range=(x_max-x_min)/8;
y_range=(y_max-y_min)/8;
if(x_range>y_range) %»√ÕºœÒ±‰∑Ω
    y_range=5*x_range-4*y_range;
else
    x_range=5*y_range-4*x_range;
end

x_min=x_min-x_range;x_max=x_max+x_range;
y_min=y_min-y_range;y_max=y_max+y_range;

%2.initialize parameters
n=length(time);

Rzs_ob=cell(n,nz);
for i=1:n
    for j=1:nz
        Rzs_ob{i,j}=Rz.*(rand(3)+0.4);%random error for covariance
    end
end

Rzs_ob=rzs;
Rzs=Rzs_ob;      %º∆À„∑Ω≤Ó
H=[eye(3),zeros(3,7)];

As=cell(n,1);
Qs=cell(n,1);

%%  AdaTE
sitar.Da=0.0001;
sitar.Dt=0.0005*eye(3,3);
sitar.Dt(2,2)=0.0001;
sitar.Dt(3,3)=0.0001;
sitar.alpha=0.01;
sitar.beta=0.2;

if strcmp(trajType,'snake')
    sitar.Da=0.0002;
    sitar.Dt=0.00005*eye(3,3);
end

ada_traj=zeros(10,n);       %πÏº£
ada_error=zeros(1,n);
ada_time = zeros(1,n);

dX=zeros(10,n);
cut_t=0;

flg=ones(n,1);
if strcmp(trajType,'snake')
    flg(250:300)=0;
end

for i=1:n
    % estimation
    tic;
    [As(1:i),Qs(1:i),Rzs(1:i,:),cut_t,ada_traj(:,1:i),dX(:,1:i),preX]=AdaTE(As(1:i-1),Qs(1:i-1),Rzs_ob(1:i,:),Rzs(1:i,:),H,cut_t,ada_traj(:,1:i-1),dX(:,1:i-1),obs(:,1:i,:),sitar,time);
    ada_error(i)=sqrt(mean(sum((ada_traj(1:3,1:i)-real(1:3,1:i)).^2,1)));
    
    ada_time(i) = toc;
    
    %     drawTrajectory(OAO_traj(1:3,1:i),1,x_min,x_max,y_min,y_max); %πÏº£Õº
    %     drawObserve(obs(1:3,1:i,:),2,x_min,x_max,y_min,y_max);%π€≤‚Õº
    
    %   drawObserveV2(obs(1:3,1:i,:),flg,12,x_min,x_max,y_min,y_max,i);%π€≤‚Õº
    drawObserveSimple(obs(1:3,1:i,:),flg,12,x_min,x_max,y_min,y_max);%π€≤‚Õº
%     drawTrajectory(OAO_traj(1:3,1:i),11,x_min,x_max,y_min,y_max); %πÏº£Õº
    drawTrajectorySimple(ada_traj(1:3,1:i),11,x_min,x_max,y_min,y_max); %πÏº£Õº
end
% draw result
drawObserveV2(obs(1:3,1:n,:),flg,12,x_min,x_max,y_min,y_max,0);%π€≤‚Õº
xlabel('');ylabel('');
set(gcf,'Position',[100,200,360,300]);
drawTrajectory(ada_traj(1:3,1:n),11,x_min,x_max,y_min,y_max); %πÏº£Õº
xlabel('');ylabel('');
set(gcf,'Position',[400,200,360,300]);

figure(95);
plot(1:n,ada_time*1000);
xlabel('Time ID');
ylabel('Time consumption (ms)');
title('Time consumption during online estimation');
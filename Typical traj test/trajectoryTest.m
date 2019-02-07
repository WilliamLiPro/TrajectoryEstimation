%experiment 1£º
%cruise trajectory experiment
%input:
% trajectory type: trajType
% observation number at same time: nz
% observation covariance: Rz
% colored noise range: crn_st,crn_ft
% rate of colored noise: ratio
%output:
function result=trajectoryTest(trajType,nz,Rz,crn_st,crn_ft,ratio)
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
if(x_range>y_range) %ÈÃÍ¼Ïñ±ä·½
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

if strcmp(trajType,'snake')
    Rzs_ob=rzs;
end
Rzs=Rzs_ob;      %¼ÆËã·½²î
H=[eye(3),zeros(3,7)];

As=cell(n,1);
Qs=cell(n,1);

%%  OAO
sitar.Da=0.0001;
sitar.Dt=0.00005*eye(3,3);
sitar.Dt(2,2)=0.00005;
sitar.Dt(3,3)=0.00005;
sitar.alpha=0.01;
sitar.beta=0.4;

if strcmp(trajType,'snake')
    sitar.Da=0.0002;
    sitar.Dt=0.0001*eye(3,3);
end

OAO_traj=zeros(10,n);       %¹ì¼£
OAO_error=zeros(1,n);

dX=zeros(10,n);
cut_t=0;

flg=ones(n,1);
if strcmp(trajType,'snake')
    flg(250:300)=0;
end

for i=1:n
    % ¹À¼Æ
    [As(1:i),Qs(1:i),Rzs(1:i,:),cut_t,OAO_traj(:,1:i),dX(:,1:i),preX]=OAOestimation(As(1:i-1),Qs(1:i-1),Rzs_ob(1:i,:),Rzs(1:i,:),H,cut_t,OAO_traj(:,1:i-1),dX(:,1:i-1),obs(:,1:i,:),sitar,time);
    OAO_error(i)=sqrt(mean(sum((OAO_traj(1:3,1:i)-real(1:3,1:i)).^2,1)));
    
%     drawTrajectory(OAO_traj(1:3,1:i),1,x_min,x_max,y_min,y_max); %¹ì¼£Í¼
%     drawObserve(obs(1:3,1:i,:),2,x_min,x_max,y_min,y_max);%¹Û²âÍ¼

% if i>=399
%     drawObserveV2(obs(1:3,1:i,:),flg,12,x_min,x_max,y_min,y_max);%¹Û²âÍ¼
%     xlabel('');ylabel('');
%     drawTrajectory(OAO_traj(1:3,1:i),11,x_min,x_max,y_min,y_max); %¹ì¼£Í¼
%     xlabel('');ylabel('');
% end
end

drawObserveV2(obs(1:3,1:i,:),flg,12,x_min,x_max,y_min,y_max);%¹Û²âÍ¼
xlabel('');ylabel('');
set(gcf,'Position',[100,200,360,300]);
drawTrajectory(OAO_traj(1:3,1:i),11,x_min,x_max,y_min,y_max); %¹ì¼£Í¼
xlabel('');ylabel('');
set(gcf,'Position',[400,200,360,300]);

%%  MAP
sitar2.Da=0.0000001;
sitar2.Dt=0.0000001*eye(3,3);

MAP_traj=zeros(14,n);       %¹ì¼£
MAP_error=zeros(1,n);

for i=1:n
    % ¹À¼Æ
    MAP_traj(:,1:i)=MAPestimation(MAP_traj(:,1:i),obs(:,1:i,:),Rzs_ob(1:i,:),sitar2,time);
    MAP_error(i)=sqrt(mean(sum((MAP_traj(1:3,1:i)-real(1:3,1:i)).^2,1)));
    
%     drawTrajectory(MAP_traj(1:3,1:i),1,x_min,x_max,y_min,y_max); %¹ì¼£Í¼
%     drawObserve(obs(1:3,1:i,:),2,x_min,x_max,y_min,y_max);%¹Û²âÍ¼
end
drawTrajectory(MAP_traj(1:3,1:i),1,x_min,x_max,y_min,y_max); %¹ì¼£Í¼

%%  EKF-RTS Acc
Qs{1}=diag([4,4,4,1,1,1,0.5,0.5,0.5]);
Qt=diag([4,4,4,1,1,1,0.5,0.5,0.5])/500;

palpha=0.01;
pbeta=0.1;

EKF_acc=zeros(9,n);       %¹ì¼£
EKF_rt=zeros(9,n);       %¹ì¼£
EKF_acc_error=zeros(1,n);

for i=1:n
    cut_t=max(1,i-50);
    
    % ¹À¼Æ
    if(i==1)
        EKF_acc(1:3,1)=mean(obs(:,1,:),3);
        EKF_rt(:,1)=EKF_acc(:,1);
    else
        [EKF_acc(:,1:i),Qs(1:i)]=EKF_RTS(obs(:,1:i,:),H(:,1:9),Rzs_ob,Qs(1:i-1),Qt,EKF_acc(:,1:i-1),EKF_rt(:,1:i-1),palpha,pbeta,time,cut_t,'Acc');
    end
    
    EKF_rt(:,i)=EKF_acc(:,i);
    EKF_acc_error(i)=sqrt(mean(sum((EKF_acc(1:3,1:i)-real(1:3,1:i)).^2,1)));
    
%     drawTrajectory(EKF_acc1(1:3,1:i),1,x_min,x_max,y_min,y_max); %¹ì¼£Í¼
%     drawObserve(obs(1:3,1:i,:),2,x_min,x_max,y_min,y_max);%¹Û²âÍ¼
end
drawTrajectory(EKF_acc(1:3,1:i),1,x_min,x_max,y_min,y_max); %¹ì¼£Í¼

%%  EKF-RTS DPT
Qs{1}=diag([4,4,4,1,1,1,0.4,0.2,0.2,0.2]);
Qt=diag([4,4,4,1,1,1,0.4,0.2,0.2,0.2])/500;

EKF_dpt=zeros(10,n);       %¹ì¼£
EKF_rt=zeros(10,n);       %¹ì¼£
EKF_dpt_error=zeros(1,n);

for i=1:n
    cut_t=max(1,i-50);
    
    % ¹À¼Æ
    if(i==1)
        EKF_dpt(:,1)=0.05*(rand(10,1)-0.5);
        EKF_dpt(7,1)=abs(EKF_dpt(7,1))*10;
        EKF_dpt(1:3,1)=mean(obs(:,1,:),3);
        EKF_rt(:,1)=abs(EKF_dpt(:,1));
    else
        [EKF_dpt(:,1:i),Qs(1:i)]=EKF_RTS(obs(:,1:i,:),H,Rzs_ob,Qs(1:i-1),Qt,EKF_dpt(:,1:i-1),EKF_rt(:,1:i-1),palpha,pbeta,time,cut_t,'DPT');
    end
    
    EKF_rt(:,i)=EKF_dpt(:,i);
    cur_error=sqrt(mean(sum((EKF_dpt(1:3,1:i)-real(1:3,1:i)).^2,1)));
    if isnan(cur_error)||isinf(cur_error)
        cur_error=10000;
    end
    EKF_dpt_error(i)=cur_error;
    
%     drawTrajectory(EKF_dpt1(1:3,1:i),1,x_min,x_max,y_min,y_max); %¹ì¼£Í¼
%     drawObserve(obs(1:3,1:i,:),2,x_min,x_max,y_min,y_max);%¹Û²âÍ¼
end
drawTrajectory(EKF_dpt(1:3,1:i),1,x_min,x_max,y_min,y_max); %¹ì¼£Í¼

%%  UKF-URTS with Acc
Qs{1}=diag([4,4,4,1,1,1,0.5,0.5,0.5]);
Qt=diag([4,4,4,1,1,1,0.5,0.5,0.5])/500;

UKF_acc=zeros(9,n);       %UKF-RTS with Acc
UKF_acc_error=zeros(1,n);  %UKF-RTS with Acc RMSE

UKF_rt=zeros(9,n);         %real_time localization

for i=1:n
    cut_t=max(1,i-50);
    
    % ÂË²¨
    tic;
    if(i==1)
        UKF_acc(1:3,1)=mean(obs(:,1,:),3);
        UKF_rt(:,1)=UKF_acc(:,1);
    else
        [UKF_acc(:,1:i),Qs(1:i)]=UKF_URTS(obs(:,1:i,:),H(:,1:9),Rzs_ob,Qs(1:i-1),Qt,UKF_acc(:,1:i-1),UKF_rt(:,1:i-1),palpha,pbeta,time,cut_t,'Acc');
    end
    
    UKF_rt(:,i)=UKF_acc(:,i);   %ÊµÊ±¹À¼Æ½á¹û
    
    % Í³¼ÆÎó²î
    UKF_acc_error(i)=sqrt(mean(sum((UKF_acc(1:3,1:i)-real(1:3,1:i)).^2,1)));
    
%     drawTrajectory(UKF_acc1(1:3,1:i),1,x_min,x_max,y_min,y_max); %¹ì¼£Í¼
%     drawObserve(obs(1:3,1:i,:),2,x_min,x_max,y_min,y_max);%¹Û²âÍ¼
end
drawTrajectory(UKF_acc(1:3,1:i),1,x_min,x_max,y_min,y_max); %¹ì¼£Í¼

%%  UKF-URTS with DPT
Qs{1}=diag([4,4,4,1,1,1,0.4,0.2,0.2,0.2]);
Qt=diag([4,4,4,1,1,1,0.4,0.2,0.2,0.2])/500;

UKF_dpt=zeros(10,n);       %UKF-RTS with PTM
UKF_dpt_error=zeros(1,n);   %UKF-RTS with PTM RMSE

UKF_rt=zeros(10,n);

for i=1:n
    cut_t=max(1,i-50);
    
    % ÂË²¨
    tic;
    if(i==1)
        UKF_dpt(:,1)=0.05*(rand(10,1)-0.5);
        UKF_dpt(7,1)=abs(UKF_dpt(7,1))*10;
        UKF_dpt(1:3,1)=mean(obs(:,1,:),3);
        UKF_rt(:,1)=UKF_dpt(:,1);
    else
        [UKF_dpt(:,1:i),Qs(1:i)]=UKF_URTS(obs(:,1:i,:),H,Rzs_ob,Qs(1:i-1),Qt,UKF_dpt(:,1:i-1),UKF_rt(:,1:i-1),palpha,pbeta,time,cut_t,'DPT');
    end
    UKF_rt(:,i)=UKF_dpt(:,i);   %ÊµÊ±¹À¼Æ½á¹û
    
    % Í³¼ÆÎó²î
    UKF_dpt_error(i)=sqrt(mean(sum((UKF_dpt(1:3,1:i)-real(1:3,1:i)).^2,1)));
    
%     drawTrajectory(UKF_dpt1(1:3,1:i),1,x_min,x_max,y_min,y_max); %¹ì¼£Í¼
%     drawObserve(obs(1:3,1:i,:),2,x_min,x_max,y_min,y_max);%¹Û²âÍ¼
end
drawTrajectory(UKF_dpt(1:3,1:i),1,x_min,x_max,y_min,y_max); %¹ì¼£Í¼

disp(['cruise trajectory has finished']);

%% statistics of RMSE
% observation RMSE
result.obRMSE=sqrt(mean(sum((obs(1:3,1:i)-real(1:3,1:i)).^2,1)));

% normalization
OAO_error=OAO_error/result.obRMSE;
MAP_error=MAP_error/result.obRMSE;
EKF_acc_error=EKF_acc_error/result.obRMSE;
EKF_dpt_error=EKF_dpt_error/result.obRMSE;
EKF_dpt_error(EKF_dpt_error>1)=1;
UKF_acc_error=UKF_acc_error/result.obRMSE;
UKF_dpt_error=UKF_dpt_error/result.obRMSE;
UKF_dpt_error(UKF_dpt_error>1)=1;

result.error=zeros(6,3);
result.RMSE=cell(6);

% OAO
result.RMSE{1}=OAO_error;
result.error(1,1)=mean(OAO_error); %mean
result.error(1,2)=min(OAO_error); %min
result.error(1,3)=max(OAO_error);%max

% MAP
result.RMSE{2}=MAP_error;
result.error(2,1)=mean(MAP_error); %mean
result.error(2,2)=min(MAP_error); %min
result.error(2,3)=max(MAP_error);%max

% EKF-Acc
result.RMSE{3}=EKF_acc_error;
result.error(3,1)=mean(EKF_acc_error); %mean
result.error(3,2)=min(EKF_acc_error); %min
result.error(3,3)=max(EKF_acc_error);%max

% EKF-DPT
result.RMSE{4}=EKF_dpt_error;
result.error(4,1)=mean(EKF_dpt_error); %mean
result.error(4,2)=min(EKF_dpt_error); %min
result.error(4,3)=max(EKF_dpt_error);%max

% UKF-Acc
result.RMSE{5}=UKF_acc_error;
result.error(5,1)=mean(UKF_acc_error); %mean
result.error(5,2)=min(UKF_acc_error); %min
result.error(5,3)=max(UKF_acc_error);%max

% UKF-DPT
result.RMSE{6}=UKF_dpt_error;
result.error(6,1)=mean(UKF_dpt_error); %mean
result.error(6,2)=min(UKF_dpt_error); %min
result.error(6,3)=max(UKF_dpt_error);%max
end
%path planning test
% path planning in random obstacle environment
%input:
% key points and arrive time
% real time obstacles
% key points covariance: Rzs_kp
% barriers covariance: Rzs_br
%output:
% online path plan

addpath('path planning');

%% 1.create key points and obstacles
%(1) key points
nkp=10;
time=1:10:nkp*10-9;
Zkp=zeros(3,nkp);
Rzs_kp=cell(nkp,1);

for i=1:nkp
    Zkp(:,i)=[i*60/nkp-(i-nkp/2)^2/5+rand(1)*50;50+30*sin(i/2)+rand(1)*40-20;(i-nkp/2)^2/5+30+rand(1)*20];
    Rzs_kp{i}=diag([1,1,0.25]);
end

%(1) barriers
nbr=4;%number of barriers
Zbr=zeros(3,nbr);
Rzs_br=cell(nbr,1);

for i=1:nbr
    Zbr(:,i)=(Zkp(:,2*i)+Zkp(:,2*i+1))/2+rand(3,1)*2-1;
    
    Rzs_br{i}=(rand(3)-0.5+diag([1,1,1]))*4;
    Rzs_br{i}=Rzs_br{i}+Rzs_br{i}';
end

H=[eye(3),zeros(3,6)];

% Rzs_br=cell(0);
% Zbr=[];

n=time(nkp);
As=cell(n,1);
Qs=cell(n,1);

%%  AdaTE
sitar.Da=0.2;
sitar.Dt=0.5*eye(3,3);
sitar.Dt(2,2)=0.5;
sitar.Dt(3,3)=0.1;
sitar.alpha=0.01;
sitar.beta=0.4;

ada_traj=zeros(9,n);       %trajectory
ada_error=zeros(1,n);

dX=zeros(9,n);
cut_t=0;

%Initialize
nx_a=time(2);
[As(1:nx_a),Qs(1:nx_a),cut_t,ada_traj(:,1:nx_a),dX(:,1:nx_a),preX]=AdaTE_PathPlan(As(1:nx_a),Qs(1:nx_a),Rzs_kp,Rzs_br,H,cut_t,ada_traj(:,1:nx_a),dX(:,1:nx_a),Zkp,Zbr,sitar,time,1);
% drawPathPlan(OAO_traj(:,1:nx_a),cut_t,Rzs_kp,Rzs_br,Zkp,Zbr,1,0,100,0,100) 
    
cut_t=cut_t+1;
for i=1:nkp-2
    if i<=nkp-2
        nx_a=time(i+2);
    else
        nx_a=time(i+1);
    end
    
    for j=time(i):time(i+1)-1
        [As(1:nx_a),Qs(1:nx_a),cut_t,ada_traj(:,1:nx_a),dX(:,1:nx_a),preX]=AdaTE_PathPlan(As(1:nx_a),Qs(1:nx_a),Rzs_kp,Rzs_br,H,cut_t,ada_traj(:,1:nx_a),dX(:,1:nx_a),Zkp,Zbr,sitar,time,1);
%         Zbr(:,j,:);
         drawPathPlan(ada_traj(:,1:nx_a),cut_t,Rzs_kp,Rzs_br,Zkp,Zbr,1,0,100,0,100) 
        cut_t=cut_t+1;
    end
end

xlabel('');ylabel('');
set(gcf,'Position',[100,200,360,300]);
drawPathPlan(ada_traj(:,1:nx_a),cut_t,Rzs_kp,Rzs_br,Zkp,Zbr,21,0,100,0,100) %path planned
xlabel('');ylabel('');
set(gcf,'Position',[400,200,360,300]);
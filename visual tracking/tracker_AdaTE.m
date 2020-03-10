%%  Copyright all Reserved by the Author and the Institution
%%  Author: William Li
%%  Email:  WilliamLi_Pro@163.com

function result=tracker_AdaTE(tracker,dataset,camPara,priDist,drawresult)
%   Visual tracker with the method of adaptive trajectory estimation
% . The observation for tracking is the result of fDSST
%
%   Input:
%   - tracker (structure): A visual tracker with multi-observations
%   - dataset (cell or structure): Sequence images of dataset and the ground truth 
%   - camPara (structure): Parameters of camPara, include (cx, cy, fx, fy)
%   - priDist (array): The primary distance of object, usually set 100
%
%   Output
%   - result (structure): The result of experiment, include the result of:
%   online filter, trajectory and their RMSE/

if(isempty(drawresult)) %默认不绘图
    drawresult=0;
end

if(isempty(tracker))
    error('No tracker available ') ;
end
if(isempty(dataset))
    error('No dataset available ') ;
end

% 读入数据集与轨迹真值
images=cell(dataset.imageNumber,1);
for i=1:dataset.imageNumber
    images{i}=imread(dataset.imagePath{i});
end

groundTruth=cell(dataset.imageNumber,1);

for i=1:dataset.imageNumber
    cur_R=dataset.groundTruth{i};   %将区域范围变为 (x_min,x_max,y_min,y_max) 格式
    switch(dataset.gtType)
        case 'VOT'
            ys=sort([cur_R(2),cur_R(4),cur_R(6),cur_R(8)]);
            xs=sort([cur_R(1),cur_R(3),cur_R(5),cur_R(7)]);
            groundTruth{i}=[mean(xs(1:2)),mean(xs(3:4)),mean(ys(1:2)),mean(ys(3:4))];
        case 'XYrange'
            groundTruth{i}=[cur_R(1),cur_R(2),cur_R(3),cur_R(4)];
        case 'YXrange'
            groundTruth{i}=[cur_R(3),cur_R(4),cur_R(1),cur_R(2)];
        case 'points xy'
            x_min=min([cur_R(1),cur_R(3)]);
            x_max=max([cur_R(1),cur_R(3)]);
            y_min=min([cur_R(2),cur_R(4)]);
            y_max=max([cur_R(2),cur_R(4)]);
            groundTruth{i}=[x_min,x_max,y_min,y_max];
        case 'points yx'
            x_min=min([cur_R(2),cur_R(4)]);
            x_max=max([cur_R(2),cur_R(4)]);
            y_min=min([cur_R(1),cur_R(3)]);
            y_max=max([cur_R(1),cur_R(3)]);
            groundTruth{i}=[x_min,x_max,y_min,y_max];
        case 'x-w-y-h'  % 左上角x-宽度-左上角y-高度
            groundTruth{i}=[cur_R(1),cur_R(1)+cur_R(2),cur_R(3),cur_R(3)+cur_R(4)];
        case 'x-y-w-h'  % 左上角x-左上角y-宽度-高度
            groundTruth{i}=[cur_R(1),cur_R(1)+cur_R(3),cur_R(2),cur_R(2)+cur_R(4)];
        otherwise
            error('Form of ground truth is not available ') ;
    end
end

% 滤波参数初值
n=dataset.imageNumber;  %时刻数
time=[1:n+1];

As=cell(n,1);
Qs=cell(n,1);

Rc=diag([1,1,0.25]);	%当前方差

H=[eye(3),zeros(3,6)];  %观测矩阵

sitar.Da=10*eye(3);
% sitar.Dt=0.5*eye(3,3);
% sitar.Dt(2,2)=0.5;
% sitar.Dt(3,3)=0.1;
sitar.alpha=0.01;
sitar.beta=0.2;
cut_t=1;

oao_traj=zeros(9,n);       %轨迹
oao_traj_plane=zeros(n,4);  %平面轨迹

mtraj_overlap=ones(1,n);    %轨迹重叠率
oao_traj_RMSE=zeros(1,n);   %轨迹RMSE

dX=zeros(9,n);

%%  第一帧图像
% 跟踪器初始化
pa=groundTruth{1};  % 位置初值
width=(pa(2)-pa(1))*priDist/camPara.fx; %目标宽度
rate=(pa(4)-pa(3))/(pa(2)-pa(1));       %高宽比
pos=[(pa(3)+pa(4))/2,(pa(1)+pa(2))/2];  %中心坐标
pc=planeToCamCood(camPara,pa,width);    %坐标转换：相平面投影->相机坐标系

oao_traj_plane(1,:)=pa;     %平面位置初值
oao_traj(1:3,1)=pc(1:3)';oao_traj(4:6,1)=[0;0;0];  %轨迹初值

dx=sqrt(Rc(1,1));dy=sqrt(Rc(2,2));dz=sqrt(Rc(3,3));
seach_r=camCoodtoPlane(camPara,[pc(1:3),dx*4+1],dy/dx); %采样范围共2倍标准差
seach_r=[seach_r,(pc(3)-dz)/pc(3),(pc(3)+dz)/pc(3)];    %加上尺寸变化

tracker.frame=1;
tracker=tracker.update(tracker,images{1},pos);  %跟踪器初始化
[tracker,pos]=tracker.estimate(tracker,images{i},pos);

% 观测值及方差
ob_size=size(tracker.observe,1);    %观测个数
obs=zeros(3,n,ob_size);
Rzs_ob=cell(n,ob_size);     %采样方差(观测先验)
Rzs=cell(n,ob_size);        %采样方差(滤波校正)

for j=1:ob_size
    cam_ob=planeToCamCood(camPara,tracker.observe(j,:),width);
    obs(:,1,j)=pc(1:3);
end
[~,~,Rzs_ob(1,:)]=observeCov(obs(:,1,:),tracker.obScore,Rc);
Rzs(1,:)=Rzs_ob(1,:);

i=1;    %状态初始化
[As(1:i),Qs(1:i),Rzs(1:i,:),cut_t,oao_traj(:,1:i),dX(:,1:i),~]=AdaTEAcc(As(1:i-1),Qs(1:i-1),Rzs_ob(1:i,:),Rzs(1:i,:),H,cut_t,oao_traj(:,1:i-1),dX(:,1:i-1),obs(:,1:i,:),sitar,time);

%%  后续更新
for i=2:n
    tracker.frame=i;
    
    % 预测
    dt=time(i+1)-time(i);
%     pre_x=AMLpredictACC(mlf_traj(:,i-1),alpha,dt);
%     pa=camCoodtoPlane(camPara,[pre_x(1:3)',width],rate);
%     pos=[(pa(3)+pa(4))/2,(pa(1)+pa(2))/2];  %中心坐标
    
    % 跟踪观测
%     dx=sqrt(Rc(1,1));dy=sqrt(Rc(2,2));dz=sqrt(Rc(3,3));
%     seach_r(1:4)=camCoodtoPlane(camPara,[pre_x(1:3)',dx*4+1],dy/dx);  %采样范围共4倍标准差
%     seach_r(5:6)=[(pre_x(3)-dz)/pre_x(3),(pre_x(3)+dz)/pre_x(3)];  %加上尺寸变化
    [tracker,pos]=tracker.estimate(tracker,images{i},pos);  %跟踪器观测更新
    
    n_ob=size(tracker.observe,1);
    for j=1:n_ob
        cam_ob=planeToCamCood(camPara,tracker.observe(j,:),width);
        obs(:,i,j)=cam_ob(1:3)';
    end
    [~,~,Rzs_ob(i,:)]=observeCov(obs(:,i,:),tracker.obScore,Rc);
    Rzs(i,:)=Rzs_ob(i,:);
    
    % 滤波校正
    [As(1:i),Qs(1:i),Rzs(1:i,:),cut_t,oao_traj(:,1:i),dX(:,1:i),~]=AdaTEAcc(As(1:i-1),Qs(1:i-1),Rzs_ob(1:i,:),Rzs(1:i,:),H,cut_t,oao_traj(:,1:i-1),dX(:,1:i-1),obs(:,1:i,:),sitar,time);
    
    % 坐标转换
    for j=1:i
        oao_traj_plane(j,:)=camCoodtoPlane(camPara,[oao_traj(1:3,j)',width],rate);
    end
    tracker=tracker.update(tracker,images{i},pos);
%     disp(pos);
    
    % 统计中心误差
    sum_e_traj=0;
    for j=1:i
        e_traj=oao_traj_plane(j,1:4)-groundTruth{j};
        ex=mean(e_traj(1:2));
        ey=mean(e_traj(3:4));
        sum_e_traj=sum_e_traj+ex^2+ey^2;
    end

    oao_traj_RMSE(i)=sqrt(sum_e_traj/i);
    
    % 统计重叠率
    sum_overlap_traj=0; %轨迹重叠率之和
    for j=1:i
        cur_overlap=overlap(oao_traj_plane(j,1:4),groundTruth{j});
        sum_overlap_traj=sum_overlap_traj+cur_overlap;
    end
    mtraj_overlap(i)=sum_overlap_traj/i;    %轨迹平均重叠率
    
%     sum_overlap_rt=sum_overlap_rt+cur_overlap;
%     mrt_overlap(i)=sum_overlap_rt/i;        %实时平均重叠率
    
    % 绘图
    [im_rows,im_cols,~]=size(images{i});
    if(drawresult)
        drawTrackResult(images{i},oao_traj(:,1:i),oao_traj_plane(i,:),camPara,1,0,im_cols,0,im_rows,0);
    end
    
    if(mod(i,10)==0)
        disp(['Current image: ',int2str(i)]);
    end
end

%   输出结果
result.trajectory=oao_traj_plane;       %结果
result.TrajectoryError=oao_traj_RMSE;   %中心误差RMSE
result.TrajectoryOverlap=mtraj_overlap; %平均重叠率
end


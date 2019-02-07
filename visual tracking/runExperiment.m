%%  Copyright all Reserved by Producer
%    Author: William Li
%    Email:  WilliamLi_Pro@163.com

%   测试在线 Online Adaptive Optimization算法
%   对比算法：EKF-RTS with Acc,EKF-RTS with PTM,UKF-URTS with Acc,UKF-URTS with PTM
%   数据：VOT数据

%   输出结果：
%   观测图 轨迹图 实时滤波图 RMSE-时间 overlap-时间
%   RMSE统计表(各个滤波算法的实时滤波误差RMSE，与轨迹RMSE的均值与方差)
%   Overlap统计表

%%  add path
cd('E:\物体检测与识别\论文工作\基于多模态融合的目标跟踪\程序\实验20180926\视觉跟踪');
addpath('../基本函数','../../相关滤波跟踪');

%%  define tracker and dataset

%   初始化数据集
im_path='E:\Dataset\Tracking Dataset\VOT2013\备份\Singer1\img';
gt_path='E:\Dataset\Tracking Dataset\VOT2013\备份\Singer1';
gt_type='x-y-w-h';  %数据集 Ground truth 格式

dataset=loadDataset(im_path,gt_path,gt_type);

%   摄像机参数
image0=imread(dataset.imagePath{1});
[im_n,im_m,~]=size(image0);
camPara.cy=im_n/2;
camPara.cx=im_m/2;
camPara.fy=min(im_n,im_m);
camPara.fx=camPara.fy;

%   目标距离
priDist=100;

p0=dataset.groundTruth{1};

%%  参数设置

% initialize the tracker
params.padding = 2.0;                   % extra area surrounding the target
params.output_sigma_factor = 1/16;		% standard deviation for the desired translation filter output
params.scale_sigma_factor = 1/16;       % standard deviation for the desired scale filter output
params.lambda = 1e-2;					% regularization weight (denoted "lambda" in the paper)
params.interp_factor = 0.025;			% tracking model learning rate (denoted "eta" in the paper)
params.num_compressed_dim = 18;         % the dimensionality of the compressed features
params.refinement_iterations = 1;       % number of iterations used to refine the resulting position in a frame
params.translation_model_max_area = inf;% maximum area of the translation model
params.interpolate_response = 1;        % interpolation method for the translation scores
params.resize_factor = 1;               % initial resize

params.number_of_scales = 17;           % number of scale levels
params.number_of_interp_scales = 33;    % number of scale levels after interpolation
params.scale_model_factor = 1.0;        % relative size of the scale sample
params.scale_step = 1.02;               % Scale increment factor (denoted "a" in the paper)
params.scale_model_max_area = 512;      % the maximum size of scale examples
params.s_num_compressed_dim = 'MAX';    % number of compressed scale feature dimensions

pos = [p0(1),p0(2)];
target_sz = [p0(3),p0(4)];
params.init_pos = floor(pos([2 1])) + floor(target_sz([2 1])/2);
params.wsize = floor(target_sz([2 1]));
params.s_frames = dataset.imagePath;

%% 参数生成
[tracker,pos]=paramsInitialize(params);

%   定义跟踪器
tracker.estimate=@DSST_estimate;
tracker.update=@DSST_update;
tracker.sampleNumber=1;

%%  OAO-DPT
oao_result=tracker_OAO(tracker,dataset,camPara,priDist,0);

%%  pure trcker
pure_result=trackerRun(tracker,dataset,camPara,priDist,0);

%%  draw the overlap
frame_label=[1:dataset.imageNumber];

figure(2);
plot(frame_label,pure_result.TrajectoryOverlap,'m:',frame_label,oao_result.TrajectoryOverlap,'c');

y_min=min([pure_result.TrajectoryOverlap,oao_result.TrajectoryOverlap]);
y_max=max([pure_result.TrajectoryOverlap,oao_result.TrajectoryOverlap]);

axis([1,dataset.imageNumber,y_min,y_max]);
title('VOT-Singer1');
set(gcf,'Position',[500,200,300,250]);

re_compare=[pure_result.TrajectoryOverlap(dataset.imageNumber),oao_result.TrajectoryOverlap(dataset.imageNumber)];

%%  draw the box comparing
frameid=3;
figure(3);
imshow(imread(dataset.imagePath{frameid}));

%GT
rectangle('Position',dataset.groundTruth{frameid}, 'EdgeColor','y');   %gt
%DSST
box_p=pure_result.trajectory(frameid,:);
rect_position_vis = [box_p([1,3]) , box_p([2,4])-box_p([1,3])+1];
rectangle('Position',rect_position_vis, 'EdgeColor','b');   %gt

%DSST+OAO
box_p=oao_result.trajectory(frameid,:); %real time
rect_position_vis = [box_p([1,3]) , box_p([2,4])-box_p([1,3])+1];
rectangle('Position',rect_position_vis, 'EdgeColor','r');   %gt

%% save image
set(gca,'position',[0 0 1 1]);
print(gcf,'-djpeg','-loose',['截图/Walking ',num2str(frameid),'.jpg']);
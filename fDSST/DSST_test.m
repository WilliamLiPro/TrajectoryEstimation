%%  变尺寸相关滤波测试函数
%   DSST test

%%  读取数据
im_path='E:\Dataset\Tracking Dataset\VOT2013\备份\Tiger1\img';
gt_path='E:\Dataset\Tracking Dataset\VOT2013\备份\Tiger1';
gt_type='x-y-w-h';
dataset=loadDataset(im_path,gt_path,gt_type);

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
tracker.sampleNumber=7;

%%  开始实验
for i=1:dataset.imageNumber
    % 读取图像
    image=imread(dataset.imagePath{i});
    tracker.frame=i;
%     if(i==2)
%         image=imread(dataset.imagePath{1});
%     end
    
    % 估计位置
    if(i>1)    %第一幅图像
        [tracker, pos] = DSST_estimate( tracker ,image ,pos );
    end
    
    % 更新参数
    tracker = DSST_update( tracker ,image ,pos);
    
    %显示图像
    rect_position_vis = [pos([2,1]) - tracker.target_sz([2,1])/2, tracker.target_sz([2,1])];
    if i == 1
        figure(1);
        im_handle = imshow(image, 'Border','tight', 'InitialMag', 100 + 100 * (length(image) < 500));
        rect_handle = rectangle('Position',rect_position_vis, 'EdgeColor','g');
        text_handle = text(10, 10, int2str(i));
        set(text_handle, 'color', [0 1 1]);
    else
        try
            set(im_handle, 'CData', image)
            set(rect_handle, 'Position', rect_position_vis)
            set(text_handle, 'string', int2str(i));
            
        catch
            return
        end
    end
    drawnow;
end
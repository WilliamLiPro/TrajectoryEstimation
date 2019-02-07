function [tracker,pos]=paramsInitialize(basicPara)
% fDSST 参数初始化函数
% 输入：
%   @basicPara 基本参数 
% 输出：
%   @tracker 初始化跟踪器
%   @pos 初始化目标坐标

padding = basicPara.padding;   %初始化
output_sigma_factor = basicPara.output_sigma_factor;   %初始化
lambda = basicPara.lambda; %estimate
interp_factor = basicPara.interp_factor;   %train
refinement_iterations = basicPara.refinement_iterations;   %estimate
translation_model_max_area = basicPara.translation_model_max_area; %初始化
nScales = basicPara.number_of_scales;  %初始化, estimate, train
nScalesInterp = basicPara.number_of_interp_scales; %初始化，estimate
scale_step = basicPara.scale_step; %初始化
scale_sigma_factor = basicPara.scale_sigma_factor;%初始化
scale_model_factor = basicPara.scale_model_factor;%初始化
scale_model_max_area = basicPara.scale_model_max_area;%初始化
interpolate_response = basicPara.interpolate_response;%estimate
num_compressed_dim = basicPara.num_compressed_dim;%初始化，train

tracker.lambda=lambda;
tracker.interp_factor=interp_factor;
tracker.refinement_iterations=refinement_iterations;
tracker.nScales=nScales;
tracker.nScalesInterp=nScalesInterp;
tracker.interpolate_response=interpolate_response;
tracker.num_compressed_dim=num_compressed_dim;


s_frames = basicPara.s_frames;%存储图像
pos = floor(basicPara.init_pos);%目标位置，各个都用
target_sz = floor(basicPara.wsize * basicPara.resize_factor); %当前目标尺寸，为时变量

tracker.s_frames=s_frames;
tracker.target_sz=target_sz;

num_frames = numel(s_frames);%图像个数
tracker.num_frames=num_frames;

init_target_sz = target_sz;%初始目标尺寸，为定量

if prod(init_target_sz) > translation_model_max_area
    currentScaleFactor = sqrt(prod(init_target_sz) / translation_model_max_area);
else
    currentScaleFactor = 1.0;
end
tracker.currentScaleFactor=currentScaleFactor;

% target size at the initial scale
base_target_sz = target_sz / currentScaleFactor;%基本目标尺寸，为定量
tracker.base_target_sz=base_target_sz;

%window size, taking padding into account
sz = floor( base_target_sz * (1 + padding ));   %搜索框初始尺寸，定量
tracker.sz=sz;

featureRatio = 4;   %特征尺寸，定量
tracker.featureRatio=featureRatio;

output_sigma = sqrt(prod(floor(base_target_sz/featureRatio))) * output_sigma_factor;%定量
use_sz = floor(sz/featureRatio);    %压缩后搜索框尺寸，仅用于初始化
rg = circshift(-floor((use_sz(1)-1)/2):ceil((use_sz(1)-1)/2), [0 -floor((use_sz(1)-1)/2)]);
cg = circshift(-floor((use_sz(2)-1)/2):ceil((use_sz(2)-1)/2), [0 -floor((use_sz(2)-1)/2)]);

[rs, cs] = ndgrid( rg,cg);
y = exp(-0.5 * (((rs.^2 + cs.^2) / output_sigma^2)));
yf = single(fft2(y));%该段代码仅用于计算目标卷积结果（高斯分布）

tracker.y=y;
tracker.yf=yf;

interp_sz = size(y) * featureRatio; %estimate用
tracker.interp_sz=interp_sz;

cos_window = single(hann(floor(sz(1)/featureRatio))*hann(floor(sz(2)/featureRatio))' );%该段代码仅用于计算cos遮罩
tracker.cos_window=cos_window;

if nScales > 0
    scale_sigma = nScalesInterp * scale_sigma_factor;
    
    scale_exp = (-floor((nScales-1)/2):ceil((nScales-1)/2)) * nScalesInterp/nScales;
    scale_exp_shift = circshift(scale_exp, [0 -floor((nScales-1)/2)]);
    
    interp_scale_exp = -floor((nScalesInterp-1)/2):ceil((nScalesInterp-1)/2);
    interp_scale_exp_shift = circshift(interp_scale_exp, [0 -floor((nScalesInterp-1)/2)]);
    
    scaleSizeFactors = scale_step .^ scale_exp;%用于 estimate,train
    interpScaleFactors = scale_step .^ interp_scale_exp_shift;%用于train
    
    tracker.scaleSizeFactors=scaleSizeFactors;
    tracker.interpScaleFactors=interpScaleFactors;
    
    ys = exp(-0.5 * (scale_exp_shift.^2) /scale_sigma^2);
    ysf = single(fft(ys));%尺寸变化train用的目标高斯分布
    scale_window = single(hann(size(ysf,2)))';%尺寸变化 estimate,train 用的cos遮罩
    
    tracker.ysf=ysf;
    tracker.scale_window=scale_window;
    
    %make sure the scale model is not to large, to save computation time
    if scale_model_factor^2 * prod(init_target_sz) > scale_model_max_area
        scale_model_factor = sqrt(scale_model_max_area/prod(init_target_sz));
    end
    
    %set the scale model size
    scale_model_sz = floor(init_target_sz * scale_model_factor);%尺寸变化 estimate,train 用的模型尺寸
    tracker.scale_model_sz=scale_model_sz;
    
    im = imread(s_frames{1});
    
    %force reasonable scale changes 尺寸参数范围estimate用
    min_scale_factor = scale_step ^ ceil(log(max(5 ./ sz)) / log(scale_step));%estimate
    max_scale_factor = scale_step ^ floor(log(min([size(im,1) size(im,2)] ./ base_target_sz)) / log(scale_step));%estimate
    
    tracker.min_scale_factor=min_scale_factor;
    tracker.max_scale_factor=max_scale_factor;
    
    max_scale_dim = strcmp(basicPara.s_num_compressed_dim,'MAX');  %最大尺寸维度 train 用
    if max_scale_dim
        s_num_compressed_dim = length(scaleSizeFactors);    % train 用
    else
        s_num_compressed_dim = basicPara.s_num_compressed_dim;
    end
    tracker.max_scale_dim=max_scale_dim;
    tracker.s_num_compressed_dim=s_num_compressed_dim;
end

% initialize the projection matrix
projection_matrix = []; %estimate,train 使用
tracker.projection_matrix=projection_matrix;

rect_position = zeros(num_frames, 4);%目标框范围，estimate输出结果
tracker.rect_position=rect_position;

end
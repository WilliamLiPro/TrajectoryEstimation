function [tracker, pos_new] = DSST_estimate( tracker ,image ,pos )
%DSST_estimate: Discriminative Scale Space Tracking 论文算法复现 - 位置与尺度估计部分
%   输入：
%   - tracker 目标跟踪器
%   - image 当前图像
%   - pos 目标上一时刻位置

%   输出：
%   - tracker 目标跟踪器
%   - pos_new 目标当前位置

% 数据准备
refinement_iterations = tracker.refinement_iterations;
sz = tracker.sz;
currentScaleFactor = tracker.currentScaleFactor;
projection_matrix = tracker.projection_matrix;
cos_window = tracker.cos_window;
featureRatio = tracker.featureRatio;   %特征缩放尺度
interp_sz = tracker.interp_sz;

hf_num = tracker.hf_num;    %训练结果：分子
hf_den = tracker.hf_den;    %训练结果：分母
lambda = tracker.lambda;    %松弛参数

interpolate_response = tracker.interpolate_response;

old_pos = inf(size(pos));
iter = 1;

%translation search
while iter <= refinement_iterations && any(old_pos ~= pos)
    [xt_npca, xt_pca] = get_subwindow(image, pos, sz, currentScaleFactor);
    
    xt = feature_projection(xt_npca, xt_pca, projection_matrix, cos_window);
    xtf = fft2(xt);
    
    responsef = sum(hf_num .* xtf, 3) ./ (hf_den + lambda);
    
    % if we undersampled features, we want to interpolate the
    % response so it has the same size as the image patch
    if interpolate_response > 0
        if interpolate_response == 2
            % use dynamic interp size
            interp_sz = floor(size(y) * featureRatio * currentScaleFactor);
        end
        
        responsef = resizeDFT2(responsef, interp_sz);
    end
    
    response = ifft2(responsef, 'symmetric');
    
    [row, col] = find(response == max(response(:)), 1);
    disp_row = mod(row - 1 + floor((interp_sz(1)-1)/2), interp_sz(1)) - floor((interp_sz(1)-1)/2);
    disp_col = mod(col - 1 + floor((interp_sz(2)-1)/2), interp_sz(2)) - floor((interp_sz(2)-1)/2);
    
    switch interpolate_response
        case 0
            translation_vec = round([disp_row, disp_col] * featureRatio * currentScaleFactor);
        case 1
            translation_vec = round([disp_row, disp_col] * currentScaleFactor);
        case 2
            translation_vec = [disp_row, disp_col];
    end
    
    old_pos = pos;
    pos = pos + translation_vec;
    
    iter = iter + 1;
end
pos_new = pos;

%scale search
% 参数初始化
base_target_sz = tracker.base_target_sz;
scaleSizeFactors = tracker.scaleSizeFactors;
scale_model_sz = tracker.scale_model_sz;
scale_basis = tracker.scale_basis;
scale_window = tracker.scale_window;
sf_num = tracker.sf_num;
sf_den = tracker.sf_den;
nScalesInterp = tracker.nScalesInterp;
nScales=tracker.nScales;

%尺度搜索
if nScales > 0
    
    %create a new feature projection matrix
    [xs_pca, xs_npca] = get_scale_subwindow(image,pos,base_target_sz,currentScaleFactor*scaleSizeFactors,scale_model_sz);
    
    xs = feature_projection_scale(xs_npca,xs_pca,scale_basis,scale_window);
    xsf = fft(xs,[],2);
    
    scale_responsef = sum(sf_num .* xsf, 1) ./ (sf_den + lambda);
    
    interp_scale_response = ifft( resizeDFT(scale_responsef, nScalesInterp), 'symmetric');
    
    recovered_scale_index = find(interp_scale_response == max(interp_scale_response(:)), 1);
    
    %set the scale
    currentScaleFactor = currentScaleFactor * tracker.interpScaleFactors(recovered_scale_index);
    %adjust to make sure we are not to large or to small
    min_scale_factor=tracker.min_scale_factor;
    max_scale_factor=tracker.max_scale_factor;
    if currentScaleFactor < min_scale_factor
        currentScaleFactor = min_scale_factor;
    elseif currentScaleFactor > max_scale_factor
        currentScaleFactor = max_scale_factor;
    end
    
    tracker.currentScaleFactor=currentScaleFactor;  %尺度更新
end

% 更新观测量
basic_rangey=[-floor((tracker.sampleNumber-1)/2):floor((tracker.sampleNumber-1)/2)];
basic_rangex=[-floor((tracker.sampleNumber-1)/2):floor((tracker.sampleNumber-1)/2)];

rangey=mod(basic_rangey+disp_row, interp_sz(1))+1;
rangex=mod(basic_rangex+disp_col, interp_sz(2))+1;

tracker.obScore=reshape(response(rangey,rangex),[tracker.sampleNumber^2,1]);
[oby,obx]=meshgrid(basic_rangey*currentScaleFactor,basic_rangex*currentScaleFactor);
tracker.observe=zeros(tracker.sampleNumber^2,4);
tracker.observe(:,1)=reshape(obx,[tracker.sampleNumber^2,1])+pos(2);
tracker.observe(:,3)=reshape(oby,[tracker.sampleNumber^2,1])+pos(1);

target_sz=floor(base_target_sz * currentScaleFactor);

cur_sz=target_sz(2:-1:1);
tracker.observe(:,1)=tracker.observe(:,1)-floor((cur_sz(ones(tracker.sampleNumber^2,1),1)-1)/2);
tracker.observe(:,2)=tracker.observe(:,1)+cur_sz(ones(tracker.sampleNumber^2,1),1)-1;
tracker.observe(:,3)=tracker.observe(:,3)-floor((cur_sz(ones(tracker.sampleNumber^2,1),2)-1)/2);
tracker.observe(:,4)=tracker.observe(:,3)+cur_sz(ones(tracker.sampleNumber^2,1),2)-1;
end
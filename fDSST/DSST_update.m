function tracker = DSST_update( tracker ,image ,pos )
%DSST_update: Discriminative Scale Space Tracking 论文算法复现 - 参数更新部分
%   输入：
%   - tracker 目标跟踪器
%   - image 当前图像
%   - pos 目标上一时刻位置
%   - GaussKernel 高斯核
%   - scaleGauss 尺度高斯核 

%   输出：
%   - tracker 目标跟踪器

% 参数初始化
sz=tracker.sz;
currentScaleFactor=tracker.currentScaleFactor;
num_compressed_dim=tracker.num_compressed_dim;
interp_factor=tracker.interp_factor;
cos_window=tracker.cos_window;
frame=tracker.frame;
%this is the training code used to update/initialize the tracker

%Compute coefficients for the tranlsation filter
[xl_npca, xl_pca] = get_subwindow(image, pos, sz, currentScaleFactor);

if frame == 1
    h_num_pca = xl_pca;
    h_num_npca = xl_npca;
    
    % set number of compressed dimensions to maximum if too many
    num_compressed_dim = min(num_compressed_dim, size(xl_pca, 2));
else
    h_num_pca = tracker.h_num_pca;
    h_num_npca = tracker.h_num_npca;
    
    h_num_pca = (1 - interp_factor) * h_num_pca + interp_factor * xl_pca;
    h_num_npca = (1 - interp_factor) * h_num_npca + interp_factor * xl_npca;
end;
tracker.h_num_pca = h_num_pca;
tracker.h_num_npca = h_num_npca;

data_matrix = h_num_pca;

[pca_basis, ~, ~] = svd(data_matrix' * data_matrix);
projection_matrix = pca_basis(:, 1:num_compressed_dim);
tracker.projection_matrix=projection_matrix;

hf_proj = fft2(feature_projection(h_num_npca, h_num_pca, projection_matrix, cos_window));
hf_num = bsxfun(@times, tracker.yf, conj(hf_proj));
tracker.hf_num = hf_num;

xlf = fft2(feature_projection(xl_npca, xl_pca, projection_matrix, cos_window));
new_hf_den = sum(xlf .* conj(xlf), 3);

if frame == 1
    hf_den = new_hf_den;
else
    hf_den = tracker.hf_den;
    hf_den = (1 - interp_factor) * hf_den + interp_factor * new_hf_den;
end
tracker.hf_den = hf_den;

%Compute coefficents for the scale filter
nScales=tracker.nScales;
base_target_sz=tracker.base_target_sz;
scaleSizeFactors=tracker.scaleSizeFactors;
scale_model_sz=tracker.scale_model_sz;
s_num_compressed_dim=tracker.s_num_compressed_dim;
scale_window=tracker.scale_window;

if nScales > 0
    
    %create a new feature projection matrix
    [xs_pca, xs_npca] = get_scale_subwindow(image, pos, base_target_sz, currentScaleFactor*scaleSizeFactors, scale_model_sz);
    
    if frame == 1
        s_num = xs_pca;
    else
        s_num = tracker.s_num;
        s_num = (1 - interp_factor) * s_num + interp_factor * xs_pca;
    end;
    tracker.s_num = s_num;
    
    bigY = s_num;
    bigY_den = xs_pca;
    
    if tracker.max_scale_dim
        [scale_basis, ~] = qr(bigY, 0);
        [scale_basis_den, ~] = qr(bigY_den, 0);
    else
        [U,~,~] = svd(bigY,'econ');
        scale_basis = U(:,1:s_num_compressed_dim);
    end
    scale_basis = scale_basis';
    tracker.scale_basis = scale_basis;
    
    %create the filter update coefficients
    sf_proj = fft(feature_projection_scale([],s_num,scale_basis,scale_window),[],2);
    sf_num = bsxfun(@times,tracker.ysf,conj(sf_proj));
    tracker.sf_num = sf_num;
    
    xs = feature_projection_scale(xs_npca,xs_pca,scale_basis_den',scale_window);
    xsf = fft(xs,[],2);
    new_sf_den = sum(xsf .* conj(xsf),1);
    
    if frame == 1
        sf_den = new_sf_den;
    else
        sf_den = tracker.sf_den;
        sf_den = (1 - interp_factor) * sf_den + interp_factor * new_sf_den;
    end;
    tracker.sf_den=sf_den;
end

target_sz = floor(base_target_sz * currentScaleFactor);
tracker.target_sz = target_sz;
end

function ft_scales=scaleFeaturs(img, pr, modelSize, scales, max_scale, cos_window, featureRatio)
%scaleFeaturs: 图像多尺度特征提取函数 提取的特征包括：RGB/gray与HOG特征
%   输入：
%   - img 当前图像
%   - pr 采样位置 (x,y)
%   - modelSize 模型尺寸
%   - scales 采样相对规模
%   - max_scale 采样最大规模
%   - cos_window cos权重遮罩
%   - featureRatio 特征采样间隔

%   输出：
%   - ft_scales 多层局部特征

n_sc=length(scales);    %采样规模数
n_feature=fix(modelSize(1)/featureRatio)*fix(modelSize(2)/featureRatio)*31; %每一层特征维数
% ft_scales_0=zeros(n_feature,n_sc);                                %出初始特征
% ft_scales=zeros(n_sc,n_sc);                                       %压缩特征
ft_scales=single(zeros(n_feature,n_sc));                            %压缩特征
% ft_size=zeros(n_sc,2);

for i=1:n_sc
    % 计算搜索范围
    cur_scale=min(scales(i),max_scale);
    
    % 该尺度下的特征
    cur_ft=regionFeatures(img, pr, cur_scale, max_scale, modelSize, cos_window, featureRatio);
    ft_scales(:,i)=reshape(cur_ft(:,:,1:31), n_feature, 1);
end

% 主成分分析(待定)
% ft_scales=
end

function [ft,range]=regionFeatures(img, pr, cur_scale, max_scale, modelSize, cos_window, featureRatio)
%regionFeatures: 图像局部特征提取函数 提取的特征包括：RGB/gray与HOG特征
%   输入：
%   - img 当前图像
%   - pr 当前位置
%   - cur_scale 当前图像比例
%   - max_scale 最大比例
%   - modelSize 模型尺寸
%   - cos_window cos权重遮罩
%   - featureRatio 特征采样间隔

%   输出：
%   - ft 局部特征
%   - sz 采样区域尺寸

% 数据准备
if(cur_scale>max_scale) %约束相对尺寸
    cur_scale=max_scale;
end

xr=floor(modelSize(1)*cur_scale); %搜索区域尺寸
yr=floor(modelSize(2)*cur_scale);

xr=max(xr,2);
yr=max(yr,2);

range=zeros(1,4);   %像素采样范围：2倍目标尺寸
range(1)=floor(pr(1)-(xr-1)/2);
range(2)=range(1)+xr-1;
range(3)=floor(pr(2)-(yr-1)/2);
range(4)=range(3)+yr-1;

% 离散采样得到尺寸归一化区域图像
xs=floor(range(1):range(2));
ys=floor(range(3):range(4));

% 限制采样范围
[n,m,~]=size(img);
ys(ys<1)=1;
ys(ys>n)=n;
xs(xs<1)=1;
xs(xs>m)=m;

% 提取局部图像并缩放
c=size(img,3);
% sub_img=imresize(img(ys,xs,:),[modelSize(2),modelSize(1)],'bilinear');
sub_img=mexResize(img(ys,xs,:),[modelSize(2),modelSize(1)], 'auto');

% 加上cos遮罩并降采样
cos_img=single(sub_img).*cos_window(:,:,ones(c,1));

% 计算特征
ft=featureExtract(cos_img,modelSize,featureRatio);

end

function ft=featureExtract(img,modelSize,featureRatio)
%featureExtract: 特征提取函数 提取的特征包括：RGB/gray与HOG特征
%   输入：
%   - img 图像
%   - modelSize 模型尺寸
%   - featureRatio 特征采样间隔

%   输出：
%   - ft 当前图像的特征图

c=size(img,3);
new_sz=fix([modelSize(2)/featureRatio,modelSize(1)/featureRatio]);
ft=single(zeros(new_sz(1),new_sz(2),31+c));

% 获取Hog特征
ft(:,:,1:32)=fhog(single(img),featureRatio);
% sp_img=imresize(img,fix([modelSize(2)/featureRatio,modelSize(1)/featureRatio]),'bilinear');
sp_img=mexResize(img,new_sz, 'auto');

% 将RGB作为特征
if(size(img,3)==3)
    ft(:,:,32:34)=single(sp_img);
else
    ft(:,:,32)=single(sp_img);
end
end
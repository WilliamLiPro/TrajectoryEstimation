%   数据集加载函数
function dataset=loadDataset(im_path,gt_path,gt_type)
%   输入：
%   - 数据集图像路径 im_path
%   - 标注路径 gt_path
%   - 图像真值存储格式 gt_type

%   输出：
%   - 数据集参数 dataset 包括：
%       - 图像路径 dataset.imagePath 
%       - 每一帧的真值 dataset.groundTruth
%       - 真值存储格式 dataset.gtType
%       - 图像个数 dataset.imageNumber

%% 读取图片
%   获取图片格式
im_type=['/*.jpg';'/*.png';'/*.bmp'];

for i=1:3
    img_path_list = dir([im_path,im_type(i,:)]); %获取该文件夹中所有jpg格式的图像
    img_num = length(img_path_list);        %获取图像总数量
    
    if(img_num)
        break;
    end
end

dataset.imagePath=cell(img_num,1);
for i=1:img_num
    dataset.imagePath{i}=fullfile(im_path,img_path_list(i).name);
end

dataset.imageNumber=img_num;
dataset.gtType=gt_type;

%%  读取Ground Truth
dataset.groundTruth=cell(img_num,1);

gt_type=gt_path(length(gt_path)-3:length(gt_path));
if(~strcmp(gt_type,'.txt'))
    flist = dir([gt_path,'/*.txt']); %获取该文件夹中所有txt
    gt_path=fullfile(gt_path,flist(1).name);
end

fpn = fopen (gt_path, 'r');           %打开文档 

id=0;
while (~feof(fpn) )
    id=id+1;
    dataset.groundTruth{id} = str2num(fgetl(fpn));
end
end
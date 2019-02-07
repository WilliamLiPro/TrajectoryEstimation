function [ s_inf ] = textSimpleTrans( s_Region,crx,cry, seta)
%textLinearTrans 纹理的简单线性状态转移函数
%   计算单个像素在给定纹理特征下周围像素对其状态的影响
%   状态量为像素灰度(0-1之间)
%输入:
%s_Region 参考点邻域状态
%crx,cry邻域中心坐标
%seta 线性推理参数
%输出：
%s_inf 推理点纹理状态

s_Region(cry,crx)=0;
s_inf=s_Region*seta;
end
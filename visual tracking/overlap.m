%%  Copyright all Reserved by Producer
%%  Author: William Li
%%  Email:  WilliamLi_Pro@163.com

%   计算两个方框之间的重叠率
function f=overlap(box1,box2)
%   Calculate the overlap of two boxes

%   Input:
%   - box1 (1*4 or 4*1 array): Box range (x_min, x_max, y_min, y_max)
%   - box2 (1*4 or 4*1 array): Box range (x_min, x_max, y_min, y_max)

%   Output
%   - f (1*1 array): The overlap of two boxes

s_box1=(box1(2)-box1(1))*(box1(4)-box1(3));
s_box2=(box2(2)-box2(1))*(box2(4)-box2(3));

x_rgmin=min([box1(2),box2(2)])-max([box1(1),box2(1)]);  %框交集
y_rgmin=min([box1(4),box2(4)])-max([box1(3),box2(3)]);
if x_rgmin<0
    x_rgmin=0;
end
if y_rgmin<0
    y_rgmin=0;
end
s_r=x_rgmin*y_rgmin;    %交集面积

f=s_r/(s_box1+s_box2-s_r);
end
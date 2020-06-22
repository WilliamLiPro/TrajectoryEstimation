%   绘制轨迹图
function drawTrajectorySimple(Xs,fig_id,x_min,x_max,y_min,y_max)
% 按时序的颜色
[~,n]=size(Xs);
if n<2
    %轨迹太短，无法绘制
    return;
end

figure(fig_id);hold off;

%定义8种颜色
color_l = ind2rgb(1:8,cool(8));

%绘制轨迹
cutid=[1,ceil((1:8)*n/8)];
has_draw=0;
for i=1:8
    if cutid(i)==cutid(i+1)
        %跳过长度为0的部分
        continue;
    end
    
    st_id = cutid(i);
    end_id = cutid(i+1);
    plot(Xs(1,st_id:end_id),Xs(2,st_id:end_id),'-','color',color_l(1,i,:));
    if has_draw==0
        hold on;
        has_draw = 1;
    end
end

axis([x_min,x_max,y_min,y_max]);
title('Trajectory');
% xlabel('x');
% ylabel('y');
% zlabel('z');

drawnow;
end
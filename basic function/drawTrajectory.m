%   绘制轨迹图
function drawTrajectory(Xs,fig_id,x_min,x_max,y_min,y_max)
% 按时序的颜色
[~,n]=size(Xs);

figure(fig_id);
clf;
hold on;

for i=1:n-1
    rate=2*i/n-1;
    
    color_l(1)=max(-rate,0);
    color_l(2)=1-abs(rate);
    color_l(3)=max(rate,0);
    
    plot3(Xs(1,i:i+1),Xs(2,i:i+1),Xs(3,i:i+1),'-','color',color_l);
end
hold off;
axis([x_min,x_max,y_min,y_max]);

xlabel('x');
ylabel('y');
zlabel('z');
end
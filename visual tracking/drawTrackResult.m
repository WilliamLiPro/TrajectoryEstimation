%% 绘制观测框和轨迹
%   rate:高宽比
function drawTrackResult(image,Xs,p_new,camera,fig_id,x_min,x_max,y_min,y_max,draw_traj)
% 按时序的颜色
[~,n]=size(Xs);

figure(fig_id);
clf;
hold on;
imshow(image);
hold off;

hold on;

% 绘制轨迹
if draw_traj
    for i=1:n-1
        rate=2*i/n-1;
        
        color_l(1)=max(-rate,0);
        color_l(2)=1-abs(rate);
        color_l(3)=max(rate,0);
        
        x1=Xs(1,i)/Xs(3,i)*camera.fx+camera.cx;y1=Xs(2,i)/Xs(3,i)*camera.fy+camera.cy;
        x2=Xs(1,i+1)/Xs(3,i+1)*camera.fx+camera.cx;y2=Xs(2,i+1)/Xs(3,i+1)*camera.fy+camera.cy;
        
        plot([x1,x2],[y1,y2],'-','color',color_l);
    end
end

% 绘制目标框
x1=p_new(1);
x2=p_new(2);
y1=p_new(3);
y2=p_new(4);
plot([x1,x2,x2,x1,x1],[y1,y1,y2,y2,y1],'y');

hold off;

pause(0.01);
%axis([x_min,x_max,y_min,y_max]);
%set(gca,'YDir','reverse');
end
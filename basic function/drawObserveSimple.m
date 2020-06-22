%   绘制观测图
% flag用于标记显示的内容
function drawObserveSimple(p_ob,flag,fig_id,x_min,x_max,y_min,y_max)
n=size(p_ob,2);

% 剔除为0的数据
p_ob(p_ob==0)=nan;

%剔除标记为空的数据
p_ob(:,flag(1:n)==0,:)=nan;

%图像id
figure(fig_id);hold off;

%定义8种颜色
color_l = ind2rgb(8:-1:1,winter(8));

cutid=[1,ceil((1:8)*n/8)];
has_draw=0;
for i=1:8
    if cutid(i)==cutid(i+1)
        %跳过长度为0的部分
        continue;
    end
    
    st_id = cutid(i);
    end_id = cutid(i+1);
    
    x=p_ob(1,st_id:end_id,:);
    y=p_ob(2,st_id:end_id,:);
    
    plot(x(:),y(:),'.','color',color_l(1,i,:),'MarkerSize',3);
    if has_draw==0
        hold on;
        has_draw = 1;
    end
end

axis([x_min,x_max,y_min,y_max]);
title('Observation');
% xlabel('x');
% ylabel('y');
% zlabel('z');

drawnow;
end
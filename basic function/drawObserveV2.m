%   绘制观测图
% flag用于标记显示的内容
function drawObserveV2(p_ob,flag,fig_id,x_min,x_max,y_min,y_max,UpTime)
n=size(p_ob,2);
lz=size(p_ob,3);

if UpTime==1 %更新时间
    if p_ob(1,UpTime)~=0
        rate=2*UpTime/n-1;
        color_l=[max(-rate,0),1-abs(rate),max(rate,0)];
        
        figure(fig_id);
        hold on;
        data1(:,:)=p_ob(:,UpTime,:);
        plot3(data1(1,:),data1(2,:),data1(3,:),'.','color',color_l,'MarkerSize',3,'LineWidth',1);
    end
    
    return;
end

% 剔除为0的数据
id=0;
for i=1:n
    if(p_ob(1,i)~=0)
        id=id+1;
    end
end
p_ob_draw=zeros(3,id,lz);
id=1;
for i=1:n
    if(p_ob(1,i)~=0)
        p_ob_draw(:,id,:)=p_ob(:,i,:);
        id=id+1;
    end
end

% 按时序着色
figure(fig_id);
clf;
hold on;

for i=1:n
    if(flag(i)==0)  %跳过0标记点
        continue;
    end
    
    rate=2*i/n-1;
    
    color_l(1)=max(-rate,0);
    color_l(2)=1-abs(rate);
    color_l(3)=max(rate,0);
    
    data1(:,:)=p_ob_draw(:,i,:);
    plot3(data1(1,:),data1(2,:),data1(3,:),'.','color',color_l,'MarkerSize',3,'LineWidth',1);
end
hold off;
axis([x_min,x_max,y_min,y_max]);

xlabel('x');
ylabel('y');
zlabel('z');
end
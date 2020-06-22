%   ���ƹ۲�ͼ
% flag���ڱ����ʾ������
function drawObserveSimple(p_ob,flag,fig_id,x_min,x_max,y_min,y_max)
n=size(p_ob,2);

% �޳�Ϊ0������
p_ob(p_ob==0)=nan;

%�޳����Ϊ�յ�����
p_ob(:,flag(1:n)==0,:)=nan;

%ͼ��id
figure(fig_id);hold off;

%����8����ɫ
color_l = ind2rgb(8:-1:1,winter(8));

cutid=[1,ceil((1:8)*n/8)];
has_draw=0;
for i=1:8
    if cutid(i)==cutid(i+1)
        %��������Ϊ0�Ĳ���
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
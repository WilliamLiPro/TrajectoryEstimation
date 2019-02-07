%绘制轨迹动态规划结果
function drawPathPlan(Xs,cut_t,Rzs_kp,Rzs_br,Zkp,Zbr,fig_id,x_min,x_max,y_min,y_max)
%平面
figure(fig_id);
clf;
hold on;

planTraj(Xs,cut_t,Rzs_kp,Rzs_br,Zkp,Zbr,x_min,x_max,y_min,y_max)

%   key point
nkp=length(Rzs_kp);
for i=1:nkp
    drawCircle(Zkp(:,i),Rzs_kp{i},-0.9,fig_id);
end

%   barier
nbr=length(Rzs_br);
for i=1:nbr
    drawCircle(Zbr(:,i),Rzs_br{i},0.9,fig_id);
end

shading flat;
lightangle(0,90);

%   path
[~,n]=size(Xs);
plot3(Xs(1,1:cut_t+1),Xs(2,1:cut_t+1),Xs(3,1:cut_t+1),'b-');
plot3(Xs(1,cut_t+1:n),Xs(2,cut_t+1:n),Xs(3,cut_t+1:n),'m-');

hold off;
axis([x_min,x_max,y_min,y_max,-1,80]);
end

function drawCircle(p,R,col,fig_id)
[dqx,dqy,dqz]=sphere(20);

rm=chol(R);

dpz=rm(3,3)*dqz+p(3);
dpy=rm(2,2)*dqy+rm(2,3)*dqz+p(2);
dpx=rm(1,1)*dqx+rm(1,2)*dqy+rm(1,3)*dqz+p(1);

figure(fig_id);
cl=col*ones(20,20);
% cl(4)=col;
surf(dpx,dpy,dpz,cl);

end

function planTraj(Xs,cut_t,Rzs_kp,Rzs_br,Zkp,Zbr,x_min,x_max,y_min,y_max)
%绘制热力图及平面轨迹
rx=x_max-x_min;
ry=y_max-y_min;
im=zeros(ry,rx);

%   key point
nkp=length(Rzs_kp);
for i=1:nkp
    cur_R=Rzs_kp{i};
    sx=2*sqrt(cur_R(1,1));
    sy=2*sqrt(cur_R(2,2));
    cx=Zkp(1,i);
    cy=Zkp(2,i);
    min_x=ceil(max([cx-sx,1]));
    max_x=ceil(min([cx+sx,rx]));
    min_y=ceil(max([cy-sy,1]));
    max_y=ceil(min([cy+sy,ry]));
    
    x=min_x:max_x;
    y=min_y:max_y;
    
    inv_R=inv(cur_R);
    hxx=exp(-(x-cx).^2*inv_R(1,1)/2);
    hyy=exp(-(y-cy).^2*inv_R(2,2)/2);
    hxy=exp(-(y-cy)'*(x-cx)*inv_R(1,2));
    h=(hyy'*hxx).*hxy;
    im(y,x)=im(y,x)-h;
end

%   barier
nbr=length(Rzs_br);
for i=1:nbr
    cur_R=Rzs_br{i};
    sx=2*sqrt(cur_R(1,1));
    sy=2*sqrt(cur_R(2,2));
    cx=Zbr(1,i);
    cy=Zbr(2,i);
    min_x=ceil(max([cx-sx,1]));
    max_x=ceil(min([cx+sx,rx]));
    min_y=ceil(max([cy-sy,1]));
    max_y=ceil(min([cy+sy,ry]));
    
    x=min_x:max_x;
    y=min_y:max_y;
    
    inv_R=inv(cur_R);
    hxx=exp(-(x-cx).^2*inv_R(1,1)/2);
    hyy=exp(-(y-cy).^2*inv_R(2,2)/2);
    hxy=exp(-(y-cy)'*(x-cx)*inv_R(1,2));
    h=(hyy'*hxx).*hxy;
    im(y,x)=im(y,x)+h;
end

%   path
[~,n]=size(Xs);
plot(Xs(1,1:cut_t+1),Xs(2,1:cut_t+1),'b-');
plot(Xs(1,cut_t+1:n),Xs(2,cut_t+1:n),'m-');

surf(im);
colormap(jet);
end
%%  坐标转换：方框的相机坐标系->相平面投影
%   采用一阶模型描述相机，不考虑畸变补偿
%   输入：
%   相机基本参数 camera (包括相机中心坐标 camera.cx camera.cy 相机焦距 camera.fx camera.fy)
%   方框位置及范围 pc=[x,y,d,w]
%   方框高宽比 rate

%   输出：
%   相平面位置 pp=[u1,u2,v1,v2]

function pp=camCoodtoPlane(camera,pc,rate)
pp=zeros(1,4);

d=pc(3);
half_w=pc(4)/2;
half_h=half_w*rate;

pp(1)=(pc(1)-half_w)/d*camera.fx+camera.cx;
pp(2)=(pc(1)+half_w)/d*camera.fx+camera.cx;
pp(3)=(pc(2)-half_h)/d*camera.fy+camera.cy;
pp(4)=(pc(2)+half_h)/d*camera.fy+camera.cy;

end
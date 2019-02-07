%%  坐标转换：方框的相平面投影->相机坐标系
%   采用一阶模型描述相机，不考虑畸变补偿
%   输入：
%   相机基本参数 camera (包括相机中心坐标 camera.cx camera.cy 相机焦距 camera.fx camera.fy)
%   方框的相平面位置 pp=[u1,u2,v1,v2]
%   方框宽度 w

%   输出：
%   相机坐标系位置 pc=[x,y,w,h]

function pc=planeToCamCood(camera,pp,w)
pc=zeros(1,4);

d=abs(w*camera.fx/(pp(2)-pp(1)));

pc(1)=((pp(1)+pp(2))/2-camera.cx)*d/camera.fx;  %x
pc(2)=((pp(3)+pp(4))/2-camera.cy)*d/camera.fy;  %y
pc(3)=d;                                        %d
pc(4)=w;                                        %w

end
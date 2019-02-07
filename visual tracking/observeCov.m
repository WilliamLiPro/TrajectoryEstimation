%%  方差生成函数
%   根据各个观测值的匹配置信度计算位置均值、总体方差与局部方差
%   输入：
%   观测位置 obs
%   匹配置信度 belief

%   输出：
%   期望位置 mp
%   总体方差 Rc
%   各个观测方差 Rzs
function [mp,Rc,Rzs]=observeCov(obs,belief,Rc)
% 计算均值
[m_ob,n_ob]=size(obs);
mp=zeros(m_ob,1);
s_blf=sum(belief);
% belief=belief/sum(belief);  %权值归一化

for i=1:n_ob
    mp=mp+obs(:,i)*belief(i);
end
mp=mp/s_blf;
% 
% % 方差
% Rc=3*Rc+0.01*eye(3);    %防止奇异
% for i=1:n_ob
%     dob=obs(:,i)-mp;
%     Rc=Rc+(belief(i)*dob)*dob';
% end
% thre=sum([Rc(1,1),Rc(2,2),Rc(3,3)])/9;
% for i=1:3
%     if(Rc(i,i)<thre)
%         Rc(i,i)=Rc(i,i)+thre;
%     else
%         if(Rc(i,i)>9*thre)
%             Rc(i,i)=Rc(i,i)-thre;
%         end
%     end
% end
% Rc=Rc/(s_blf*4);

% 局部方差
thre=-min([belief;0])+0.001;
Rzs=cell(1,n_ob);
for i=1:n_ob
    Rzs{i}=Rc/(belief(i)+thre); %防止分母为0
end
end
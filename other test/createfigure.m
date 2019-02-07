function createfigure(X1, YMatrix1)
%CREATEFIGURE(X1, YMATRIX1)
%  X1:  x 数据的矢量
%  YMATRIX1:  y 数据的矩阵

%  由 MATLAB 于 01-Nov-2018 15:49:11 自动生成

% 创建 figure
figure1 = figure;

% 创建 axes
axes1 = axes('Parent',figure1);
box(axes1,'on');
hold(axes1,'all');

% 使用 plot 的矩阵输入创建多行
plot1 = plot(X1,YMatrix1,'Parent',axes1);
set(plot1(1),'Color',[1 0.5 0.1],'DisplayName','OAO-DPT');
set(plot1(2),'Color',[0.1 0.6 1],'DisplayName','MAP-CT');

% 创建 xlabel
xlabel('k');

% 创建 ylabel
ylabel('ms');

% 创建 legend
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.184447605500234 0.745529507306775 0.381578947368421 0.156133828996283]);


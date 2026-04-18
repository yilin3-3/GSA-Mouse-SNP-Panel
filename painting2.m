
% 1. 读取CSV文件
% 注意：MATLAB 2014的csvread函数只能读取数值数据
% 需要先读取表头，再读取数值矩阵
file_name = '315669/hamming_distance.csv'
file_name = '512/hamming_distance.csv'
file_name = '297/hamming_distance.csv'
% 读取第一行作为表头（品系名称）
fid = fopen(file_name, 'r');
header_line = fgetl(fid);
fclose(fid);

% 使用textscan解析表头
header = textscan(header_line, '%s', 'Delimiter', ',');
strain_names = header{1}(2:end); % 去掉第一个空单元格

% 读取数值矩阵（跳过第一行）
distance_matrix = csvread(file_name, 1, 1);

% 2. 创建热图
figure('Position', [100, 100, 1200, 800]); % 设置图形窗口大小

% 创建自定义颜色映射
% 0为白色，色阶从1开始
% 使用jet颜色映射，但将第一个颜色设置为白色
cmap = jet(256); % 获取256色的jet颜色映射
cmap(1,:) = [1 1 1]; % 将第一个颜色设置为白色（对应值0）

% 绘制热图
imagesc(distance_matrix);
colormap(cmap); % 应用自定义颜色映射

% 3. 设置颜色条
c = colorbar;
caxis([0 max(distance_matrix(:))]); % 设置颜色范围，0对应白色

% 4. 添加表头标签
% 设置x轴和y轴标签
set(gca, 'XTick', 1:length(strain_names), 'XTickLabel', strain_names);
set(gca, 'YTick', 1:length(strain_names), 'YTickLabel', strain_names);

% 旋转x轴标签以避免重叠
xtickangle(45);

% 5. 设置背景颜色
set(gcf, 'Color', [1 1 1]); % 设置图形背景为白色
set(gca, 'Color', [1 1 1]); % 设置坐标轴背景为白色

% 6. 添加标题和标签
title('小鼠品系间距离矩阵热图', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('品系名称', 'FontSize', 12);
ylabel('品系名称', 'FontSize', 12);

% 7. 调整布局
grid off;
axis tight;

% 8. 保存图像
saveas(gcf, 'distance_heatmap.png');
print(gcf, '-dpng', '-r300', 'distance_heatmap_highres.png'); % 高分辨率保存

disp('热图已生成并保存为 distance_heatmap.png');

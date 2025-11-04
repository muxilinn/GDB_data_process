clc; clear;
%% 读取数据
% 定义文件夹路径
folderPath = 'C:\Users\Liluo\Desktop\轨道板数据\data_GDB\k13';
data = catch_data(folderPath, 1);

%% 数据预处理aa
% 有没有不用data.字段名的方法
% 直接将结构体转换为数组
dataArray = struct2array(data);
% 将结构体的字段名转换为数组
fieldNames = fieldnames(data);

time = dataArray(1,1).vib_data(:,1);
z_data = dataArray(1,1).vib_data(:,3);
% 绘制第一个变量的空间图
figure;
plot(time,z_data);
% 用fieldNames来设置标题和轴标签,横轴对了,但是纵轴没有对
% 纵轴是z_data
% 注意一下图窗,从数据的第一个点到最后一个点,图窗大小要适当调整
% 可以用axis函数来调整图窗大小
axis([time(1) time(end) min(z_data) max(z_data)]);
title([strrep(fieldNames{1},'_','\_') '轨道板振动图像']);
xlabel(strrep(fieldNames{1},'_','\_'));
ylabel(strrep('垂向振动数据','_','\_'));

%% 计算RMS值
% 设置窗口参数
window_size = 1000;
window_type = 'hann';  % 或 'hamming', 'rectangular'

% 手动实现滑动窗口RMS
% 使用示例
rms_z = sliding_rms_vectorized(z_data, window_size, window_type);
% rms_z = zeros(size(z_data));
% for i = window_size:length(z_data)
%     window_data = z_data(i-window_size+1:i);
    
%     % 应用窗函数（可选）
%     if strcmp(window_type, 'hann')
%         win = hann(window_size);
%     elseif strcmp(window_type, 'hamming')
%         win = hamming(window_size);
%     else
%         win = ones(window_size, 1);  % 矩形窗
%     end
    
%     window_data_windowed = window_data .* win;
%     rms_z(i) = sqrt(mean(window_data_windowed.^2));
% end

% 处理前window_size-1个点（设为NaN或使用较小窗口）
rms_z(1:window_size-1) = NaN;

% 并将RMS值用折线图表示
figure;
plot(time, rms_z);
title('垂向振动数据的RMS值');
xlabel(strrep(fieldNames{1},'_','\_'));
ylabel(strrep('RMS值','_','\_'));


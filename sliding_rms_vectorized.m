function rms_z = sliding_rms_vectorized(data, window_size, window_type)
    % 数据长度
    n = length(data);
    
    % 创建索引矩阵
    indices = (1:window_size)' + (0:n-window_size);
    
    % 提取所有窗口数据
    window_data = data(indices);
    
    % 应用窗函数
    if strcmp(window_type, 'hann')
        win = hann(window_size);
    elseif strcmp(window_type, 'hamming')
        win = hamming(window_size);
    else
        win = ones(window_size, 1);
    end
    
    % 向量化计算
    window_data_windowed = window_data .* win;
    rms_values = sqrt(mean(window_data_windowed.^2, 1));
    
    % 构建完整结果
    rms_z = [nan(1, window_size-1), rms_values];
end

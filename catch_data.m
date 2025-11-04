function [data] = catch_data(folderPath, columnIndex)
    % 提取该文件夹下所有的.mat文件
    fileNames = dir(fullfile(folderPath, '*.mat'));
    % 初始化一个空矩阵来存储所有数据
    % allData = [];
    dataStruct = struct();
    % 用一个for循环来遍历每个文件
    for i = 1:length(fileNames)
        % 读取当前文件的数据
        currentData = load(fullfile(folderPath, fileNames(i).name));
        % 提取指定的列数据
        columnData = currentData(:, columnIndex);
        % 获取基本文件名并创建有效变量名
        [~, baseName, ~] = fileparts(fileNames(i).name);
        % 创建有效的MATLAB标识符
        validName = matlab.lang.makeValidName(baseName);
        % 使用结构体存储所有数据（推荐方式）
        dataStruct.(validName) = columnData;
    end
    % 将所有数据赋值给输出变量data
    data = dataStruct;
end

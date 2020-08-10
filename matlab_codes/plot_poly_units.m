clear all; clc;

% fileList = dir('./mat_files/*.mat');
fileList = dir('../py_codes/*.mat');

for ct1 = 1:size(fileList,1)
    mat_name = [fileList(ct1).folder, '/',fileList(ct1).name];
    mat_name
    cc_coors_analysis(mat_name);
end
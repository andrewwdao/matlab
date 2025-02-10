% add folders to search path
cur_dir = pwd;
folders = {'doa','localisation'};
for ii = 1:length(folders)
    addpath(genpath(fullfile(cur_dir, folders{ii})));
end
clear; clc; close all;

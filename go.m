% add folders to search path
cur_dir = pwd;
folders = {'basics', 'utils', 'arrays', 'doa','localisation', 'data'};
for ii = 1:length(folders)
    addpath(genpath(fullfile(cur_dir, folders{ii})));
end
clear; clc; close all;

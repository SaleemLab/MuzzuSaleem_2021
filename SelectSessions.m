%% Tomaso Muzzu - UCL - 16/01/2019

function FileNames = SelectSessions
clear all
% load file (selectable by user)
DataFolder = ['X:\DATA\SUBJECTS'];
FileNames = uigetfile_n_dir(DataFolder,'Select recording file');

end
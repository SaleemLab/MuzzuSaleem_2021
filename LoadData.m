%% Tomaso Muzzu - UCL - 13/09/2019

% load data of visual perturbation project from the IBN-vision server 

function [ProjectData AM_UnitResponses AM_Param AM_Speed AM_Time AM_UOI AM_EyeTracking] = LoadData

clear all
[~,currHostname] = system('hostname');
currHostname = currHostname(1:end-1);
switch currHostname
    case 'saleem12'
        serverName = ['X:\ibn-vision\']; 
        DataFolder = [serverName 'DATA\PROJECTS\VisPerturbation'];     
    case 'Tommy-ThinkPad'
        serverName = ['C:\Users\Tommy\OneDrive\WORK\UCL\VisPerturbation']
        DataFolder = serverName;     
    otherwise
        serverName = ['X:\']; 
        DataFolder = [serverName 'DATA\PROJECTS\VisPerturbation'];     
end

FileNames = uigetfile_n_dir(DataFolder,'Select recording file');

for i=1:length(FileNames)
    load(FileNames{1,i});
end
if exist('ProjectDataCTRL','var') 
    ProjectData = ProjectDataCTRL;
elseif exist('ProjectDataALL','var')
    ProjectData = ProjectDataALL;
end
if ~exist('AM_UOI','var')
    AM_UOI = 1;
end    
if ~exist('AM_Time','var')
    AM_Time = 1;
end   
if ~exist('AM_EyeTracking','var')
    AM_EyeTracking = 1;
end  

if exist('AM_UnitResponsesALL','var')
    AM_UnitResponses = AM_UnitResponsesALL;
end  
if exist('AM_ParamALL','var')
    AM_Param = AM_ParamALL;
end  
if exist('AM_ParamALL','var')
    AM_Speed = AM_SpeedALL;
end  

end
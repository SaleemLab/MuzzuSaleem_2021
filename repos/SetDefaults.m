function DIRS = SetDefaults(varargin)
% Function to load up the default directories for the data or recording
% setup as defaults
% Optional input arguments are:
% ['animal', animal_name]: to get animal specific folders
% ['rig', rig_name]: to get rig specific folders (example calibration
% files)
%
% Note that all the folders have a filesep at the end


%
% Aman Saleem
% 10 Apr 2019

%% input parameters and flags
pnames = {'animal', 'rig'};
dflts  = {[], []};
[animalName, rigName] = ...
    internal.stats.parseArgs(pnames,dflts,varargin{:});
if ~isempty(rigName);     getRigDirs = true;   else;    getRigDirs = false;     end
if ~isempty(animalName);  getAnimalDirs = true;else;    getAnimalDirs = false;  end

%% Main program
[~, hostname] = system('hostname');
hostname(end) = [];

% Point to the correct drive on the computer.
switch hostname
    case 'saleem12'
        DIRS.serverName = ['X:\ibn-vision' filesep];
    otherwise
        DIRS.serverName = ['X:\ibn-vision' filesep];
        if ~exist(DIRS.serverName)
            DIRS.serverName = ['X:\'];
        end
        disp(['This machine is: ' hostname ', it is not listed and there might be errors']);
        disp(['Using this: ' DIRS.serverName]);
        disp('If there are errors, try adding a new hostname in the case - like saleem12 above');
end

DIRS.Data       = [DIRS.serverName 'DATA' filesep];
DIRS.Rigs       = [DIRS.Data 'RIGS' filesep];
DIRS.Subjects   = [DIRS.Data 'SUBJECTS' filesep];

DIRS.Database           = [DIRS.Data 'DATABASES' filesep];
DIRS.DatabaseRigs       = [DIRS.Database 'DB-Rigs' filesep];
DIRS.DatabaseSubjects   = [DIRS.Database 'DB-Subjects' filesep];
DIRS.DatabaseUsers      = [DIRS.Database 'DB-Users' filesep];

% Subject specific folders
if getAnimalDirs
    DIRS.AnimalDirectory    = [DIRS.Subjects animalName filesep];
    
    DIRS.Anatomy            = [DIRS.AnimalDirectory 'Anatomy' filesep];
    DIRS.BehaviourVideos    = [DIRS.AnimalDirectory 'BehaviourVideos' filesep];
    DIRS.BonVisionStimFiles = [DIRS.AnimalDirectory 'BonVisionStimFiles' filesep];
    DIRS.EyeTracking        = [DIRS.AnimalDirectory 'EyeTracking' filesep];
    DIRS.MiniscopeImaging   = [DIRS.AnimalDirectory 'MiniscopeImaging' filesep];
    DIRS.Records            = [DIRS.AnimalDirectory 'Records' filesep];
    DIRS.VRBehaviour        = [DIRS.AnimalDirectory 'VRBehaviour' filesep];
    DIRS.WidefieldImaging   = [DIRS.AnimalDirectory 'WidefieldImaging' filesep];
    DIRS.ePhys              = [DIRS.AnimalDirectory 'ePhys' filesep];
    DIRS.Processed          = [DIRS.AnimalDirectory 'Processed' filesep];
    
    %     Further subfields can be added here
    DIRS.Widefieldmaps      = [DIRS.WidefieldImaging filesep 'Widefieldmaps' filesep];
end

% Rig specific folders
if getRigDirs
    DIRS.RigDirectory   = [DIRS.Rigs rigName filesep];
    
    DIRS.CameraCalib    = [DIRS.RigDirectory 'CameraCalibration' filesep];
    DIRS.MonitorCalib   = [DIRS.RigDirectory 'MonitorCalibration' filesep];
    DIRS.RewardCalib    = [DIRS.RigDirectory 'RewardCalibration' filesep];
end
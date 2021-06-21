% Input:  a structure with any subset of the following fields (a user will be prompted
%                                           for whatever details that are not supplied)
%         .animal - a string
%         .series - a string or a number
%         .selectedExperiments - a cell array of strings or a column vector of numbers
%         .identifier - a string which is used in forming the output filename
%         .channelOrder - a string (e.g. 'CM32')
%         .car - true if Common Average Reference is to be computed
% 
% 2015-02 MO
function ztitch_openEphys(inputStruct)

SetDefaultDirs; % this script is in \\zserver\code\Spikes\
global DIRS;

if nargin < 1 || isempty(inputStruct)
  inputStruct.dummyField = [];
end;
if ~isstruct(inputStruct)
  error('ztitch_openEphys receives a single argument which must be a structure');
end;


%-------
%Determine animal name
if ~isfield(inputStruct, 'animal') || isempty(inputStruct.animal)
  tmp = inputdlg('Please provide the animal name');
  if isempty(tmp)
    error('No animal name provided');
  end;
  inputStruct.animal = tmp{1};
end;

%-------
%Get a valid identifier
if ~isfield(inputStruct, 'identifier') || isempty(inputStruct.identifier)
  tmp = inputdlg('Please provide the identifier to use in the output filenames');
  if isempty(tmp)
    error('No identifier provided');
  end;
  inputStruct.identifier = tmp{1};
elseif ~ischar(inputStruct.identifier)
  error('No valid identifier provided');
end;

%-------
%Determine the series
if ~isfield(inputStruct, 'series') || isempty(inputStruct.series)
  d=dir([DIRS.Cerebus filesep inputStruct.animal]);
  if isempty(d)
    error(['Directory ' DIRS.Cerebus filesep inputStruct.animal ' not found']);
  end;
  d=d(cell2mat({d(:).isdir})); % throw out the files, leave the subdirectories
  d=d(3:end); % throw out . and ..
  if isempty(d)
    error(['Directory ' DIRS.Cerebus filesep inputStruct.animal ' has no subdirectories']);
  end;
  d = {d.name};
  [series,tmp] = listdlg('PromptString','Please select the series',...
    'SelectionMode','single',...
    'ListString',d);
  if tmp < 1
    error('Error with selecting series...');
  end;
  inputStruct.series = d{series}; % series was an index into the list of choices, now it's a string
  clear d series;
elseif isnumeric(inputStruct.series)
  inputStruct.series = num2str(inputStruct.series);
end;

%-------
%Determine the experiments
if ~isfield(inputStruct, 'selectedExperiments') || isempty(inputStruct.selectedExperiments)
  d=dir([DIRS.Cerebus filesep inputStruct.animal filesep inputStruct.series]);
  if isempty(d)
    error(['Directory ' DIRS.Cerebus filesep inputStruct.animal filesep inputStruct.series ' not found']);
  end;
  d=d(cell2mat({d(:).isdir})); % throw out the files, leave the subdirectories
  d=d(3:end); % throw out . and ..
  if isempty(d)
    error(['Directory ' DIRS.Cerebus filesep inputStruct.animal filesep inputStruct.series ' has no subdirectories']);
  end;
  d = {d.name};
  [selectedExperiments,tmp] = listdlg('PromptString',{'Please select the', 'experiments'},...
    'InitialValue',1:length(d),...
    'ListString',d);
  if tmp < 1
    error(['Error with selecting experiments...']);
  end;
  inputStruct.selectedExperiments = {d{selectedExperiments}}; % selectedExperiments were indexes into the list of choices
  clear d selectedExperiments;
elseif isvector(inputStruct.selectedExperiments)
  tmp = {};
  for i = 1:length(inputStruct.selectedExperiments)
    tmp{i} = num2str(inputStruct.selectedExperiments(i));
  end;
  inputStruct.selectedExperiments = tmp;
  clear tmp i;
end;

% WHEN the experiment names are both single and double digit we will run
% into a problem of sorting them by name ?
assert(length(inputStruct.selectedExperiments) < 9);

%-------
% Verify that each experiment directory has the same channels
% and build a list of all the directories that will be unified into the .dat file
dirs2unify = [];
for e = 1:length(inputStruct.selectedExperiments)
  dirs2unify(end+1).name =  [DIRS.Cerebus filesep inputStruct.animal filesep inputStruct.series filesep inputStruct.selectedExperiments{e}];
  contFiles = dir([dirs2unify(end).name filesep '*.continuous']);
  if isempty(contFiles)
    error([DIRS.Cerebus filesep inputStruct.animal filesep inputStruct.series filesep inputStruct.selectedExperiments{e} filesep ' contains no .continuous file']);
  end;
  
  % rename _ADC files (external analog I/O) into _CH files
  for j=1:length(contFiles)
    assert(isempty(strfind(contFiles(j).name, 'AUX')), 'Why AUX channels are there?!')
    pos = strfind(contFiles(j).name, 'ADC');
    if isempty(pos)
      continue
    end
    assert(contFiles(j).name(pos+4) == '.', 'after channel number only ''.continuous'' suffix is expected')
    chN = str2num(contFiles(j).name(pos+3));
    if length(contFiles) < 80
      newName = [ contFiles(j).name(1:end-15) 'CH' num2str(chN+90) '.continuous'];
    else % for recordings of 128-512 channels
      newName = [ contFiles(j).name(1:end-15) 'CH' num2str(chN+900) '.continuous'];
    end
    assert(~exist([dirs2unify(end).name filesep newName], 'file'));
    movefile([dirs2unify(end).name filesep contFiles(j).name], ...
      [dirs2unify(end).name filesep newName]);
  end  
  contFiles = dir([dirs2unify(end).name filesep '*.continuous']);
  
  dirs2unify(end).ch = zeros(1, length(contFiles));
  for j=1:length(contFiles)      
    dirs2unify(end).ch(j) = getChannelDetails(contFiles(j).name);
  end    
  dirs2unify(end).ch = sort(dirs2unify(end).ch, 'ascend');
  dirs2unify(end).filePrefixName = contFiles(1).name(1:strfind(contFiles(1).name, 'CH')+1);
  clear contFiles datearray dateindex j;
end;
clear e

%------- Check that all the directories have exactly the same channels

for i = 2:length(dirs2unify)
  if length(dirs2unify(i).ch) ~= length(dirs2unify(1).ch)
    error('ztitch_openEphys: not all recordings have the same no. of channels, a condition that (for now) is not handled!');
  elseif max(abs(dirs2unify(i).ch - dirs2unify(1).ch)) > 0
    error('ztitch_openEphys: not all recordings have exactly the same set of channels?!');
  end;
end;
clear i fileH;

%------- Construct the name of the output .dat filename (which we are about to create) & check it does not exist already:
outputDir = [DIRS.multichanspikes filesep inputStruct.animal filesep inputStruct.series filesep];
if ~exist(outputDir, 'dir')
  mkdir(outputDir);
  disp(['Created directory ' outputDir]);
end;
% outputDatFilename = [outputDir filesep inputStruct.animal '_' inputStruct.series '_' inputStruct.identifier '.dat'];
outputDatFilename = [outputDir filesep inputStruct.animal '_s' inputStruct.series '_' inputStruct.identifier '.dat'];
outputMatFilename = strcat(outputDatFilename(1:end-4), '.mat');
if exist(outputDatFilename, 'file') || exist(outputMatFilename, 'file')
  error(['.dat/.mat file with the name ' outputDatFilename(1:end-4) ' already exists?!']);
end;

%------- Present the user the list of all the directories that will be stitched together
disp(['The following directories (in the order shown) are about to be stitched into ', outputDatFilename]);
disp('----------------------------------------------------------------------------------------------');
disp(char(dirs2unify(:).name));
disp('==============================================================================================');


CHANNELS_ORDER = 1:length(dirs2unify(1).ch);
ORDER_OPTIONS = {'None', 'CM32', 'CM2x2tet', 'CM16lin'};
if ~isfield(inputStruct, 'channelOrder') || isempty(inputStruct.channelOrder)
[selection,ok] = listdlg('PromptString', 'Select channel remapping option', ...
  'SelectionMode', 'single', ...
  'ListString', ORDER_OPTIONS);
if ~ok
  return;
end;
else
  for selection = 1:numel(ORDER_OPTIONS)
    if strcmp(ORDER_OPTIONS{selection}, inputStruct.channelOrder)
	  break;
	end;
  end;
end;
switch ORDER_OPTIONS{selection}
  case 'CM32'
    if length(dirs2unify(1).ch) < 32  % sanity check
      error('This cannot possibly be a CM32 recording?!');
    end;
    CHANNELS_ORDER = [3 4 14 5 13 6 12 7 11 8 10 9 15 16 2 1 32 31 17 18 24 23 25 22 26 21 27 20 28 19 29 30 33:length(dirs2unify(1).ch)];
  case 'CM16lin'
    if length(dirs2unify(1).ch) < 16  % sanity check
      error('This cannot possibly be a CM16 recording?!');
    end;
    CHANNELS_ORDER = [8 9 7 10 6 11 5 12 16 1 15 2 16 3 13 4 17:length(dirs2unify(1).ch)];
  case 'CM2x2tet'
    if length(dirs2unify(1).ch) < 16  % sanity check
      error('This cannot possibly be a CM16 recording?!');
    end;
    CHANNELS_ORDER = [1 3 5 15 7 9 11 13 6 8 10 12 2 4 14 16 17:length(dirs2unify(1).ch)];
  otherwise
    ; % channels order is not tempered with...
end;
clear selection ok

if ~isfield(inputStruct, 'car') || isempty(inputStruct.car)
car = questdlg('To apply CAR (common average referencing)?', ...
  'CAR Menu','Yes','No', 'No');
elseif inputStruct.car
  car = 'Yes';
else
  car = 'No';
end;

%------- Now we're sort of finally ready to get to the real work of writing the .dat file
outputFileH = fopen(outputDatFilename, 'w');
if outputFileH < 0
  error(['Failed to open ', outputDatFilename ' for writing']);
end;

lims = []; % this variable will contain the number of samples in each directory that goes into the output .dat
for f = 1:length(dirs2unify)
  data = [];
  for ch = 1:length(dirs2unify(f).ch) % open the files and load data
    fid = fopen([dirs2unify(f).name filesep dirs2unify(f).filePrefixName ...
      num2str(dirs2unify(f).ch(ch)) '.continuous'], 'r');
    if fid < 0
      error(['Error opening ', dirs2unify(f).filePrefixName, num2str(dirs2unify(f).ch(ch))]);
    end;
    
    fread(fid, 1024, 'char*1'); % header
    cnt = 1;
    samples = {};
    while true
      timestamp = fread(fid, 1, 'int64',0,'l');
      if isempty(timestamp)
        break;
      end;
      N = fread(fid, 1, 'uint16',0,'l');
      fread(fid, 1, 'uint16', 0, 'l'); % recordingNumber
      samples{cnt} = fread(fid, N, 'int16',0,'b');
      fread(fid, 10, 'char*1'); % recordmarker
      cnt = cnt+1;
    end;
    if numel(samples{end}) < 1024 % open-ephys crashed while writing last block...
      samples = samples(1:end-1);
    end;
    fclose(fid);
    samples=reshape(int16(cell2mat(samples)), [], 1);
    if isempty(data)
      data = int16(false(length(dirs2unify(1).ch), length(samples)));
    end;
    data(ch,:) = samples;
  end; % finished loading data from this directory
  
  lims(end+1) = int64(size(data,2)); %can be a pretty big number, not to loose precision with float representation, we make it int64
  data = data(CHANNELS_ORDER, :);
  if strcmp(car, 'Yes')
    if length(CHANNELS_ORDER) < 32
      carCh = int16(mean(data(1:16,:)));
      data(1:16,:) = bsxfun(@minus, data(1:16,:), carCh);
    else
      carCh = int16(mean(data(1:32,:)));
      data(1:32,:) = bsxfun(@minus, data(1:32,:), carCh);
    end;
    tmp=fwrite(outputFileH, [data; carCh], 'int16');
    if tmp ~= size(data,1) * size(data, 2) + numel(carCh)
      error(['Problem writing into ', outputDatFilename]);
    end;
  else
    tmp=fwrite(outputFileH, data, 'int16');
    if tmp ~= size(data,1) * size(data, 2)
      error(['Problem writing into ', outputDatFilename]);
    end;
  end;
  
end;
fclose(outputFileH);
if strcmp(car, 'Yes')
  CHANNELS_ORDER(end+1) = -999; % the CAR channel
end;
tmp = dir(outputDatFilename);
if strcmp(car, 'Yes') && tmp.bytes ~= sum(lims)*(length(dirs2unify(1).ch)+1)*2
  error('Some problem writing the .dat file: its actual size is incorrect!');
elseif ~strcmp(car, 'Yes') && tmp.bytes ~= sum(lims)*length(dirs2unify(1).ch)*2 % sanity check
  error('Some problem writing the .dat file: its actual size is incorrect!');
end;

%------- We save some parameters about the stitching into an acompanying .mat file:
save(outputMatFilename, 'lims', 'CHANNELS_ORDER', 'car', ...
  'dirs2unify', 'inputStruct');


function [number, type] = getChannelDetails(filenamestring)
assert(strcmp(filenamestring(end-10:end), '.continuous'))
pos1 = strfind(filenamestring, 'ADC');
pos2 = strfind(filenamestring, 'CH');

if isempty(pos1) && isempty(pos2)
  error('Bad filename (perhaps AUX channels present?!')
elseif isempty(pos1)
  type = 'CH';
  pos = pos2;
else
  type = 'ADC';
  pos = pos1;
end
switch type
  case 'ADC'
    assert(filenamestring(pos+5) == '.', 'after channel number only ''.continuous'' suffix is expected')
    number = num2str(filenamestring(pos+4));
  otherwise % 'CH'
    if filenamestring(pos+3) == '.' % 1 digit number
      number = str2double(filenamestring(pos+2));
    elseif filenamestring(pos+4) == '.' % 2 digit number
      number = str2double(filenamestring(pos+2:pos+3));
    elseif filenamestring(pos+5) == '.' % 3 digit number
      number = str2double(filenamestring(pos+2:pos+4));
    else
      error('after channel number only ''.continuous'' suffix is expected')
    end
end
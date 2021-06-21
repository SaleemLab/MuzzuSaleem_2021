function shortList = listSubjects(varargin)

% Loads a filtered list of subjects
% Usage: shortList = listSubjects(varargin)
%   varargin empty, gives list of all subjects on the server
% 
%   varargin: 'type': filters by subjects which have experiments of type specified
%   e.g.: listSubjects('type', 'Anatomy')
% 
%   varargin: 'user': filters by user specified
%   e.g.: listSubjects('user', 'Aman') will list all animals associated to Aman
% 
%   varargin: 'text': filters by having the text within the name of the subject
%   e.g.: listSubjects('text', 'R190') will list all subjects with R190 in the name
% 
% Additionally there is 'text2' and 'type2'.
% All combinations of the above types can be used.
%   Example usage: a = listSubjects('type', 'ePhys', 'user', 'Tomaso')
% 
% 'refresh' set to 1: searches through and updated the .csv files of the
% databases

% ___________
% Aman Saleem
% 10 Apr 2019

%% input parameters
pnames = {'type', 'user', 'text', 'type2', 'text2', 'refresh'};
dflts  = {[], [], [], [], [], 0};
[type, user, phrase, type2, phrase2, refresh_flag] = ...
    internal.stats.parseArgs(pnames,dflts,varargin{:});
%% Main program
DIRS = SetDefaults;
shortList = subjectList(refresh_flag);

if ~isempty(type);    type_active = 1;      else;    type_active = 0;   end
if ~isempty(type2);   type2_active = 1;     else;    type2_active = 0;  end
if ~isempty(user);    user_active = 1;      else;    user_active = 0;   end
if ~isempty(phrase);  phrase_active = 1;    else;    phrase_active = 0; end
if ~isempty(phrase2); phrase2_active = 1;   else;    phrase2_active = 0;end

if type_active
    typeSubjects = subjectListOfType(type, refresh_flag);
    shortList = shortList(contains(shortList, typeSubjects));
end
if type2_active
    type2Subjects = subjectListOfType(type2, refresh_flag);
    shortList = shortList(contains(shortList, type2Subjects));
end
if length(shortList)==0
    errordlg('No data left in shortlist');
end
if user_active
    userSubjects = getUserSubjects(user);
    shortList = shortList(contains(shortList, userSubjects));
end
if length(shortList)==0
    errordlg('No data left in shortlist');
end
if phrase_active
    shortList = shortList(contains(shortList, phrase));
end
if length(shortList)==0
    errordlg('No data left in shortlist');
end
if phrase2_active
    shortList = shortList(contains(shortList, phrase2));
end
if length(shortList)==0
    errordlg('No data left in shortlist');
end

    
%%%%%% Other functions
%% Subject List
    function Subjects = subjectList(refresh)
        % This just runs and updates the SubjectList DB file.
        %
        % AS 12/2/2019 - checked and added output option
        if nargin<1
            refresh = 0;
        end
        if refresh
            subs = dir(DIRS.Subjects);
            subsTable = struct2table(subs);
            Subjects = subsTable.name(3:end);
            writetable(cell2table(Subjects), [DIRS.DatabaseSubjects 'SubjectList.csv']);
        else
            tempA = readtable([DIRS.DatabaseSubjects 'SubjectList.csv']);
            Subjects = tempA.Subjects;
        end
    end

%% Subject List of Type
    function Subjects = subjectListOfType(type, refresh)
        % This is to update the lists of different type to be updated
        % Example Usage:
        %        subjectListOfType(type, refresh1_or_not0) (Good to run this first)
        %        subjects = subjectListOfType('VRBehaviour');
        %        subjects = subjectListOfType('VRBehaviour', 1);
        % AS: 12/2/2019
        
        if nargin<2
            refresh = 0;
        end
        if nargin<1
            error('Enter the type of folder to be updated as input');
        elseif ~ischar(type)
            error('Enter the type of folder to be updated, as a string');
        end
        
        subpath = [DIRS.DatabaseSubjects 'SubjectList.csv'];
        
        if refresh
            disp('NOTE: Make sure the subject list is up-to-date. This can be done by running refreshSubjectList.m');
            
            subjectList = readtable(subpath);
            nSubs = length(subjectList.Subjects);
            
            isPresent = false(1,nSubs);
            for iSub = 1:nSubs
                isPresent(iSub) = isfolder([DIRS.Subjects subjectList.Subjects{iSub} filesep type]);
            end
            Subjects = subjectList.Subjects(isPresent);
            if sum(isPresent)>0
                writetable(cell2table(Subjects), [DIRS.DatabaseSubjects  'SubjectList_' type '.csv']);
            else
                disp('No folders of this type, check spelling');
            end
        else
            temp = readtable([DIRS.DatabaseSubjects  'SubjectList_' type '.csv']);
            Subjects = temp.Subjects;
        end
    end

%% Subject List of User
    function userSubjects = getUserSubjects(userName)
        % This just loads the user's subject list.
        % AS 09/4/2019
        fileName = [DIRS.DatabaseUsers 'UserSubjectList_' userName '.csv'];
        if exist(fileName)
            s  = readtable(fileName, 'Delimiter', ' ');
            if size(s,1)>0
                userSubjects = s.Subjects;
            else
                errordlg('No Subjects for specified type & user')    
            end
        else
            errordlg('No Subject file for specified user')
        end
    end

end
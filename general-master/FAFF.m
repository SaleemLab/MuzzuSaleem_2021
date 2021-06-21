function [metaStruct, paths, names, DIRS] = FAFF(prefFile)
% FAFF: Filter Animals and Find Files/Paths
% [metaStruct, paths, names, DIRS] = FAFF(prefFile);
%
% This is basically a GUI to run the function findFiles(animal, exp_type,
% iseries, exp_list, exp_ext) with appropriate inputs and just pass on the
% outputs from the function

% Aman Saleem 11 Apr 2019

if nargin<1
    loadPrefFile = 0;
elseif ~exist(prefFile)
    loadPrefFile = 0;
    disp('No such preference file, starting with no filters')
    prefFile = [];
else
    loadPrefFile = 1;
    disp('Will load specified filter preferences')
end

if loadPrefFile
    load(prefFile);
else
    filterOpts.type  = [];
    filterOpts.user  = [];
    filterOpts.text  = [];
    filterOpts.type2 = [];
end

animal = [];
currLevel = 0;
animal_list = [];

series_list = [];
exp_list = [];
figHandles = [];
iseries = [];
exp_list = [];
exp_ext = [];
metaStruct = []; 
paths= [];
names= [];
DIRS= [];
animal_list = listSubjects('type',filterOpts.type,'user', filterOpts.user,'text', filterOpts.text, 'type2', filterOpts.type2);

createUI
set(figHandles.CE_RunButton, 'String', 'Return Animal List');
% set(figHandles.CE_RunButton, 'Enable', 'off');
% uiwait

% set(figHandles.CA_animalListBox, 'String', animal_list);
    function run_findFiles
        switch currLevel
            case 0
%                 Animal filters chosen
                DIRS = SetDefaults;
                names = listSubjects('type',filterOpts.type,'user', filterOpts.user,'text', filterOpts.text, 'type2', filterOpts.type2);
                metaStruct.DIRS = DIRS;
                metaStruct.path = DIRS.Subjects;
                metaStruct.ChosenLevel = 'Animal List';
                metaStruct.animal = [];
                metaStruct.exp_list = [];
                metaStruct.series = [];
                metaStruct.list = names;
            case 1
%                 Animal chosen
                [metaStruct, paths, names, DIRS]  = findFiles(animal);
            case 2
%                 Animal & Exp Type chosen
                [metaStruct, paths, names, DIRS]  = findFiles(animal, filterOpts.type);
            case 3
%                 Animal, Exp Type & Series chosen
                [metaStruct, paths, names, DIRS]  = findFiles(animal, filterOpts.type, iseries);
            case 4
%                   Animal, Exp Type, Series & Exp List chosen
                [metaStruct, paths, names, DIRS]  = findFiles(animal, filterOpts.type, iseries, exp_list);
            case 5
                %                   Animal, Exp Type, Series & Exp List chosen
                [metaStruct, paths, names, DIRS]  = findFiles(animal, filterOpts.type, iseries, exp_list, exp_ext);
        end
        close(figHandles.MainFig);
    end

    function animal_chosen
        currLevel = 1;
        updateSubjectInputs;
        animal = animal_list{figHandles.CA_animalListBox.Value};
        if isempty(filterOpts.type)
            msgbox('Choose an Experiment Type and then reselect an animal!')
        else 
            [~, ~, series_list]  = findFiles(animal, filterOpts.type);
            currLevel = 2;
            set(figHandles.CE_RunButton, 'String', 'Return Series List');
            set(figHandles.CS_SeriesListBox, 'String', series_list);
        end
    end

    function refresh_animal_list
        animal_list = listSubjects('type',filterOpts.type,...
            'user', filterOpts.user,'text', filterOpts.text, 'type2', filterOpts.type2, 'refresh',1);
    end

    function filters_type_chosen
        currLevel = 0;
        filterOpts.type = figHandles.CF_MainTypeTextBox.String;
        updateSubjectInputs;
        set(figHandles.CS_MainTypeText, 'String', filterOpts.type)
    end
    function filters_user_chosen
        filterOpts.user = figHandles.CF_UserTextBox.String;
        updateSubjectInputs;
    end
    function filters_text_chosen
        filterOpts.text = figHandles.CF_TextTextBox.String;
        updateSubjectInputs;
    end
    function filters_type2_chosen
        filterOpts.type2 = figHandles.CF_Type2TextBox.String;
        updateSubjectInputs;
    end
    function updateSubjectInputs
        animal_list = listSubjects('type',filterOpts.type,...
            'user', filterOpts.user,'text', filterOpts.text, 'type2', filterOpts.type2);
        set(figHandles.CA_animalListBox, 'String', animal_list);
        set(figHandles.CS_SeriesListBox, 'String', series_list);
        set(figHandles.CE_ExptListBox, 'String', exp_list);
        set(figHandles.CE_ExptListBox, 'min',0, 'max', length(exp_list));
    end
    function filters_clear
        filterOpts.type  = [];
        set(figHandles.CF_MainTypeTextBox, 'String', filterOpts.type);
        filterOpts.user  = [];
        set(figHandles.CF_UserTextBox, 'String', filterOpts.user);
        filterOpts.text  = [];
        set(figHandles.CF_TextTextBox, 'String', filterOpts.text);
        filterOpts.type2 = [];
        set(figHandles.CF_Type2TextBox, 'String', filterOpts.type2);
        set(figHandles.CS_MainTypeText, 'String', '(Of type...)')
        set(figHandles.CE_RunButton, 'String', 'Return Animal List');
        series_list = []; 
        exp_list = []; 
        updateSubjectInputs;
    end
    function filters_load
        filters_clear
        [fileName, filterFilePath] = uigetfile('*.mat', 'Choose preference file');
        fullPath = [filterFilePath fileName];
        load(fullPath);
        updateSubjectInputs;
        set(figHandles.CF_MainTypeTextBox, 'String', filterOpts.type);
        set(figHandles.CF_UserTextBox, 'String', filterOpts.user);
        set(figHandles.CF_TextTextBox, 'String', filterOpts.text);
        set(figHandles.CF_Type2TextBox, 'String', filterOpts.type2);
    end
    function filters_save
        displayText = 'Folder to save preference file';
        [fileName,filterFilePath] = uiputfile('*.mat');
        fullPath = [filterFilePath fileName];
        save(fullPath, 'filterOpts');
    end
    function series_chosen
        currLevel = 3;
        iseries = series_list{figHandles.CS_SeriesListBox.Value};
        [~, ~, exp_list]  = findFiles(animal, filterOpts.type, iseries);
        set(figHandles.CE_ExptListBox, 'String', exp_list)
        set(figHandles.CE_ExptListBox, 'min',0,'max', length(exp_list));

        set(figHandles.CE_RunButton, 'String', 'Return Experiment List');
    end
    function exp_chosen
        currLevel = 4;
        exp_list = exp_list{figHandles.CE_ExptListBox.Value};
        [~, ~, exp_list]  = findFiles(animal, filterOpts.type, iseries);
        set(figHandles.CE_RunButton, 'String', 'Return Experiment details');
    end
    function createUI
        figHandles.MainFig = figure('Name', 'Filter Animals & Find Files/Paths',...
                    'MenuBar', 'none', 'Toolbar', 'none',...
                    'NumberTitle', 'off',  'Units', 'normalized',...
                    'OuterPosition', [0.1 0.1 0.7 0.7]);
        % Overall figure, HBoxFlex: AllPanels
        figHandles.AllPanels = uix.HBoxFlex('Parent',figHandles.MainFig,'Spacing',10,'Padding',5);
        
        % B: Column 1: Filter info, VBoxFlex: Col_Filt
        figHandles.Col_Filt = uiextras.VBoxFlex('Parent',figHandles.AllPanels,'Spacing',10,'Padding',5);
        
        % A: Column 2: Animal info, VBoxFlex: Col_Animal
        figHandles.Col_Animal = uiextras.VBoxFlex('Parent',figHandles.AllPanels,'Spacing',10,'Padding',5);
            % 1. CA_title
            figHandles.CA_title = uicontrol('Style','text',...
                'Parent',figHandles.Col_Animal,...
                'String','Animal list:',...
                'fontsize',14, 'HorizontalAlignment','left');
            % 2. CA_animalListBox
            figHandles.CA_animalListBox = uicontrol('Style','listbox',...
                'Parent',figHandles.Col_Animal,...
                'fontsize',9, 'HorizontalAlignment','left',...
                'Callback', @(~,~) animal_chosen);
            set(figHandles.CA_animalListBox, 'String', animal_list)
            % 3. CA_RefreshButton
            figHandles.CA_RefreshButton = uicontrol('Style','pushbutton',...
                'Parent',figHandles.Col_Animal,...
                'background',[1 .8 1],'fontsize',14,...
                'String','Refresh','Callback', @(~,~) refresh_animal_list, 'Enable','on');
            
        
            % 4. CF_title
            figHandles.CF_title = uicontrol('Style','text',...
                'Parent',figHandles.Col_Filt,...
                'String','Filter list',...
                'fontsize',14, 'HorizontalAlignment','left');
            
            % 5. CF_MainTypeTitle
            figHandles.CF_MainTypeTitle = uicontrol('Style','text',...
                'Parent',figHandles.Col_Filt,...
                'String','Main Type',...
                'fontsize',14, 'HorizontalAlignment','left');
            % 6. CF_MainTypeTextBox
            figHandles.CF_MainTypeTextBox = uicontrol('Style','edit',...
                'Parent',figHandles.Col_Filt,...
                'fontsize',10,...
                'String',filterOpts.type, 'Callback', @(~,~) filters_type_chosen);
            % 7. CF_UserTitle
            figHandles.CF_UserTitle = uicontrol('Style','text',...
                'Parent',figHandles.Col_Filt,...
                'String','User',...
                'fontsize',14, 'HorizontalAlignment','left');
            % 8. CF_UserTextBox
            figHandles.CF_UserTextBox = uicontrol('Style','edit',...
                'Parent',figHandles.Col_Filt,...
                'fontsize',10,...
                'String',filterOpts.user, 'Callback', @(~,~) filters_user_chosen);
            % 9. CF_TextTitle
            figHandles.CF_TextTitle = uicontrol('Style','text',...
                'Parent',figHandles.Col_Filt,...
                'String','Text (in name)',...
                'fontsize',14, 'HorizontalAlignment','left');
            % 10. CF_TextTextBox
            figHandles.CF_TextTextBox = uicontrol('Style','edit',...
                'Parent',figHandles.Col_Filt,...
                'fontsize',10,...
                'String',filterOpts.text, 'Callback', @(~,~) filters_text_chosen);
            % 11. CF_Type2Title
            figHandles.CF_Type2Title = uicontrol('Style','text',...
                'Parent',figHandles.Col_Filt,...
                'String','Secondary Type',...
                'fontsize',14, 'HorizontalAlignment','left');
            % 12. CF_Type2TextBox
            figHandles.CF_Type2TextBox = uicontrol('Style','edit',...
                'Parent',figHandles.Col_Filt,...
                'fontsize',10,...
                'String',filterOpts.type2, 'Callback', @(~,~) filters_type2_chosen);
            
                % 13. CF_ClearButton
                figHandles.CF_ClearButton = uicontrol('Style','pushbutton',...
                'Parent',figHandles.Col_Filt,...
                'background',[1 .8 1],'fontsize',14,...
                'String','Clear Filters','Callback', @(~,~) filters_clear, 'Enable','on');
            
            % 13. CF_FilterButtonBox (HBox)
            figHandles.CF_FilterButtonBox = uix.HBox('Parent',figHandles.Col_Filt ,'Spacing',10,'Padding',5);
                % 15. CF_LoadButton
                figHandles.CF_LoadButton = uicontrol('Style','pushbutton',...
                'Parent',figHandles.CF_FilterButtonBox,...
                'background',[1 1 .8],'fontsize',14,...
                'String','Load','Callback', @(~,~) filters_load, 'Enable','on');
                % 16. CF_SaveButton
                figHandles.CF_SaveButton = uicontrol('Style','pushbutton',...
                'Parent',figHandles.CF_FilterButtonBox,...
                'background',[1 1 .8],'fontsize',14,...
                'String','Save','Callback', @(~,~) filters_save, 'Enable','on');
            uix.Empty('Parent', figHandles.Col_Filt);
            
        % C: Column 3: Series info, VBoxFlex: Col_Series
        figHandles.Col_Series = uiextras.VBoxFlex('Parent',figHandles.AllPanels,'Spacing',10,'Padding',5);
            % 17. CS_title
            figHandles.CS_title = uicontrol('Style','text',...
                'Parent',figHandles.Col_Series,...
                'String','Series list',...
                'fontsize',14, 'HorizontalAlignment','left');
            % 18. CS_MainTypeText
            figHandles.CS_MainTypeText = uicontrol('Style','text',...
                'Parent',figHandles.Col_Series,...
                'String','(Of type...)',...
                'fontsize',14, 'HorizontalAlignment','left');
            % 19. CS_SeriesListBox
            figHandles.CS_SeriesListBox = uicontrol('Style','listbox',...
                'Parent',figHandles.Col_Series,...
                'fontsize',9, 'HorizontalAlignment','left',...
                'Callback', @(~,~) series_chosen);
            set(figHandles.CS_SeriesListBox, 'String', series_list)
            
        % D: Column 4: Expt   info, VBoxFlex: Col_Expt
        figHandles.Col_Expt = uiextras.VBoxFlex('Parent',figHandles.AllPanels,'Spacing',10,'Padding',5);
            % 20. CE_title
            figHandles.CS_title = uicontrol('Style','text',...
                'Parent',figHandles.Col_Expt,...
                'String','Experiment list',...
                'fontsize',14, 'HorizontalAlignment','left');
            % 21. CE_ExptListBox
            figHandles.CE_ExptListBox = uicontrol('Style','listbox',...
                'Parent',figHandles.Col_Expt,...
                'fontsize',9, 'HorizontalAlignment','left',...
                'min',0,'max', 10,...
                'Callback', @(~,~) exp_chosen);
            set(figHandles.CE_ExptListBox, 'String', exp_list)
            % 22. CE_RunButton
            figHandles.CE_RunButton = uicontrol('Style','pushbutton',...
                'Parent',figHandles.Col_Expt,...
                'background','Green','fontsize',14,...
                'String','Load','Callback', @(~,~) run_findFiles, 'Enable','on');
            % Set sizes of the GUI
%             set(figHandles.CE_RunButton, 'min',0, 'max', length(exp_list));
            set(figHandles.Col_Animal,  'Sizes',[-1 -12 -2]);
            set(figHandles.Col_Filt,    'Sizes',[-1*ones(1,11) -4]);
            set(figHandles.Col_Series,  'Sizes',[-1 -1 -10]);
            set(figHandles.Col_Expt,    'Sizes',[-1 -9 -3]);
    end
end
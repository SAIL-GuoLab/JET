% A sample script for using the JET pipeline.
%
% Chen "Raphael" Liu (cl3760@columbia.edu) & Jia Guo (jg3400@columbia.edu),
% 05/14/2020.

clear
clc

%% Load the current folder and subfold
% Add all subfolders to the system path.
addpath(genpath('./'));

%% Define the subject folder.
%*************************************************************************%
% Define the path to the study data. Each study folder may contain multiple
% subject folders, under which there may be one or more scan folders (the
% multiple scan folders must come from the same subject and same study).
%*************************************************************************%
% This step can be pop-up screen. Will need to explain the hierarchy:
% Study/subject_acquisition/session.

%study_pathname = './Data/Siemens/3T/EmilyStudy/';
study_pathname = './Data/Bruker/Biospec_94_20_usr/Invivo_Mouse_Thalamus/';

%*************************************************************************%
% Search for all the subject folders of the study.
% We need to make sure that we remove '.' (the current directory) and 
% '..' (the parent directory) from the search results.
%*************************************************************************%
subject_foldername_all = dir(study_pathname);
subject_foldername_all = ...
    subject_foldername_all(...
    ~ismember({subject_foldername_all.name}, {'.', '..'}));

%% Initialization, load data, spectrum registration, and spectral fitting
% Iterate over the subjects.
for current_subject_folder_index = 1 : length(subject_foldername_all)
    JET_STR.Loader.current_subject_folder_index = current_subject_folder_index;
    
    %% Step 1. Initialization.
    if current_subject_folder_index == 1
        % A. User-defined parameters. Future option: can be loaded from an
        % excel. JET excel filename shall be built into the folder name,
        % such that if the filename is changed, the path to save the
        % results will change. 
        JET_STR.Loader.study_pathname = study_pathname;
        JET_STR.Loader.subject_foldername_all = subject_foldername_all;
        
        % B. System-defined parameters.
        JET_STR = JET_Initialization(JET_STR);
    end

    %% A. Step 2. Load data (Loader).
    JET_STR = JET_Loader(JET_STR);

    %% B. Step 3. Spectrum Registration (SR).
    disp(strcat(['>>>> Performing spectrum registration for subject: ', ...
        JET_STR.Loader.subject_foldername_all(JET_STR.Loader.current_subject_folder_index).name, '...']));

    JET_STR = JET_Spectrum_Registration(JET_STR);

    disp(strcat('Spectrum registration for subject: ', ...
        JET_STR.Loader.subject_foldername_all(JET_STR.Loader.current_subject_folder_index).name, ' >>>>>> Done!'));

    %% C. Step 4. Spectra fitting (SF).
    disp('>>>>>> Spectra fitting is running...')

    JET_STR = JET_Spectral_Fitting(JET_STR);
    
    % We need another structure, BasisSet_STR that contain the basis set
    % and fitting frequency bound.
    % In the future, shall integrate BasisSet_STR into the initialization.

    disp('...Spectra fitting is Done!')

    %% D. Step 5. Generate Report (Report).
    % Merge into all three other modules (Loader, SR, SF)!!!
    disp('>>>>>> Generating report...') 

    JET_STR = JET_Report(JET_STR);

    disp('...JET is Done!')
    
    pause(1)
end

% Remove used variables.
clear current_subject_folder_index
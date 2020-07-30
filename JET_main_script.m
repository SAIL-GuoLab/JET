% A sample script for using the JET pipeline.
%
% Chen "Raphael" Liu (cl3760@columbia.edu) & Jia Guo (jg3400@columbia.edu),
% 02/24/2020.

% 20200526 Meeting Notes
% 1. Load and process Florence KA data (Did Florence Sham instead.)
% 2. Load and process MRS_struct_GE3T_Cell (done!)
% 3. 1st and 2nd level phase matching and correction (done!)
% 4. use fitting to help addional dephase for display (done!)
% 5. SR code debug? (done!)

clear
clc

%% Load the current folder and subfold
% Add all subfolders to the system path.
addpath(genpath('./'));

%% Option 1. Use our Loader. This requires you to define the subject folder.
%*************************************************************************%
%*** Please comment out this section if you are not using this option. ***%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
JET_STR.Loader.option = 'JET';
%*************************************************************************%
% Define the path to the study data. Each study folder may contain multiple
% subject folders, under which there may be one or more scan folders (the
% multiple scan folders must come from the same subject and same study).
%*************************************************************************%
% This step can be pop-up screen. Will need to explain the hierarchy:
% scanner/Study/subject_acquisition/session.

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% !!!!!!!!!! Please change the pathname to your Study folder. !!!!!!!!!
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%study_pathname = './Data/Bruker/Biospec_94_20_usr/Invivo_Mouse_Thalamus/'; %volume coil
study_pathname = './Data/Bruker/Biospec_94_30_usr/Jia_GABA/'; %multi-array surface coil
%Define the coil type.
JET_STR.Loader.coil_type = 'multi_array_surface_coil'; % 'volume_coil'/'multi_array_surface_coil'
%Define the basisset folder
JET_STR.SF.params.basisset_foldername='MEGAPRESS_400MHz_20ppm_32768_TE68';%'MEGAPRESS_400MHz_10ppm_32768_TE68'/'MEGAPRESS_400MHz_20ppm_32768_TE68'

%% Option 2. Directly use the data loaded from Gannet Loader.
%*************************************************************************%
%*** Please comment out this section if you are not using this option. ***%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% JET_STR.Loader.option = 'Gannet';
% Gannet_filename = './Data/MRS_struct/MRS_struct_GE3T_cell.mat';
% Gannet_file = load(Gannet_filename, 'MRS_struct_GE3T_cell');
% Gannet_file = Gannet_file.MRS_struct_GE3T_cell;
% JET_STR.Loader.Gannet_file = Gannet_file;
% clear Gannet_filename
% JET_STR.Loader.coil_type = 'volume_coil'; % 'volume_coil'/'multi_array_surface_coil'

%% Initialization, load data, spectrum registration, and spectral fitting
if strcmp(JET_STR.Loader.option, 'Gannet')
    study_path_splitted = split(Gannet_file{1, 1}.gabafile{1}, '/');
    study_path_name_cell = join(study_path_splitted(1 : end - 2), '/');
    JET_STR.Loader.study_pathname = strcat(study_path_name_cell{1, 1}, '/');
    clear study_path_splitted study_path_name_cell
    for file_index = 1 : length(Gannet_file)
        study_path_splitted = split(Gannet_file{1, file_index}.gabafile{1}, '/');
        JET_STR.Loader.subject_foldername_all(file_index).name = study_path_splitted{end-1};
    end

elseif strcmp(JET_STR.Loader.option, 'JET')
    % Search for all the subject folders of the study.
    % We need to make sure that we remove '.' (the current directory) and
    % '..' (the parent directory) from the search results.
    JET_STR.Loader.subject_foldername_all = dir(study_pathname);
    JET_STR.Loader.subject_foldername_all = ...
        JET_STR.Loader.subject_foldername_all(...
        ~ismember({JET_STR.Loader.subject_foldername_all.name}, {'.', '..'}));
    JET_STR.Loader.study_pathname = study_pathname;
end

%% Iterate over the subjects.
for current_subject_folder_index = 1 : length(JET_STR.Loader.subject_foldername_all)
    JET_STR.Loader.current_subject_folder_index = current_subject_folder_index;
    %% Step 1. Initialization (Initialization).
    if current_subject_folder_index == 1
        % 1.1. User-defined parameters. Future option: can be loaded from an
        % excel. JET excel filename shall be built into the folder name,
        % such that if the filename is changed, the path to save the
        % results will change.

        % 1.2. System-defined parameters.
        JET_STR = JET_Initialization(JET_STR);
    end

    %% Step 2. Load data (Loader).    
    disp('>>>>>> Spectra loader is running...')
    JET_STR = JET_Loader(JET_STR);
    disp('>>>>>> Spectra loader is done!')
    
    %% Step 3. Spectrum Registration (SR).
    disp(strcat(['>>>> Performing spectrum registration for subject: ', ...
        JET_STR.Loader.subject_foldername_all(JET_STR.Loader.current_subject_folder_index).name, '...']));

    JET_STR = JET_Spectrum_Registration(JET_STR);

    disp(strcat('>>>>>> Spectrum registration for subject: ', ...
        JET_STR.Loader.subject_foldername_all(JET_STR.Loader.current_subject_folder_index).name, ' is done!'));

    %% Step 4. Spectra fitting (SF).
    disp('>>>>>> Spectra fitting is running...')

    JET_STR = JET_Spectral_Fitting(JET_STR); % Jia's comment: apply additional phase correction for Off, On, Diff mean spectra, s.t. report can display a phase corrected real spectra

    % We need another structure, BasisSet_STR that contain the basis set
    % and fitting frequency bound.
    % In the future, shall integrate BasisSet_STR into the initialization.

    disp('>>>>>> Spectra fitting is done!')

    %% Step 5. Generate Report (Report).
    disp('>>>>>> Generating report...')

    JET_STR = JET_Report(JET_STR);

    disp('>>>>>> JET is done!')
    
    pause(1)
end

% Remove used variables.
clear current_subject_folder_index
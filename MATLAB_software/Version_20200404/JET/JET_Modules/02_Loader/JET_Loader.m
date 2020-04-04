function JET_STR = JET_Loader(JET_STR, current_subject_folder_index)
%JET_LOADER First step in the JET pipeline: use the Loader to load data.
%
% Input arguments
% - JET_STR : The JET struct for the study.
% - current_subject_folder_index : The index of the current folder for the
% current subject. 1 = first folder, 2 = second folder, etc.
%
% Output arguments
% - JET_STR : The modified JET struct for the study.
%
% Chen "Raphael" Liu (cl3760@columbia.edu) & Jia Guo (jg3400@columbia.edu),
% 02/24/2020.

close all;

% Display the progress.
disp('>>>>>> Loading Scans...');

% Pre-initialize the loader parameters. These values will be changed
% once they are obtained from the headers of the actual files.
JET_STR = JET_Loader_preinitialize(JET_STR);

% Start loading the header files and update the parameters.
JET_STR = JET_Loader_update_params(JET_STR, current_subject_folder_index);

disp('... Scans Loaded!');
end

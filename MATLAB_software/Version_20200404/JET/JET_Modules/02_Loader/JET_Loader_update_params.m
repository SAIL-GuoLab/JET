function JET_STR = JET_Loader_update_params(JET_STR, current_subject_folder_index)
%JET_LOADER_UPDATE_PARAMS A helper function to facilitate JET_Loader, which
%is the first step of the JET pipeline. This helper function comes after 
%JET_Loader_preinitialize, and it aims to overwrite some assumed parameters
%as the true values are acquired when loading the actual file headers.
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

switch JET_STR.Loader.vendor
     case 'Bruker'
        JET_STR = JET_Loader_load_Bruker(JET_STR, current_subject_folder_index);
        
        % Here loads the Full FID Data.
        JET_STR.SR.data.organized_fid_data = JET_STR.Loader.data.organized_fid_data;
        %JET_STR.SR.data.organized_fid_data = JET_STR.Loader.data.organized_fid_data .* ...
        %    JET_STR.Loader.receiver_gain_all(current_subject_folder_index);
end

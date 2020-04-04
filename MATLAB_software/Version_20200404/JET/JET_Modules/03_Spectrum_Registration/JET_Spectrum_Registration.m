function JET_STR = JET_Spectrum_Registration(JET_STR, current_subject_folder_index)
%JET_Spectrum_Registration Second step in the JET pipeline: Perform
% spectrum registration based on the selected protocol.
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

% Create a report directory to save the figures.
JET_STR.Report.Subject_foldername = JET_STR.Loader.subject_foldername_all(current_subject_folder_index).name;
DateString = datestr(now, 'yyyy-mm-dd');
JET_STR.Report.report_dir = strcat(JET_STR.Report.report_prefix, DateString, '/', JET_STR.Report.Subject_foldername);

if (exist(JET_STR.Report.report_dir, 'dir') ~= 7)
    mkdir(JET_STR.Report.report_dir)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Apply appropriate pre-processing. %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Keep two versions of the FID data. Version 1 will be more "authentic"
%%% but with lower SNR, which will eventually be used for quantification.
%%% Version 2 will be more processed, less "authentic" but with higher SNR,
%%% which will be used as a high-SNR representation for registration.
%%%
%%% Version 1: fid_data_less_processing
%%% Only use zero-filling.
%%%
%%% Version 2: fid_data_more_processing
%%% Use zero-filling, LB, SNR enhancement, denoising.
%%%
%%% The overall goal is to calculate the transform using Version 2 and
%%% apply that transform on Version 1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Zero filling.
JET_STR.SR.data.time_zeropad = (1 : 1 : JET_STR.SR.params.ZeroFillTo) / JET_STR.SR.params.sw;

if JET_STR.SR.params.ZeroFillTo < size(JET_STR.SR.data.organized_fid_data, 1)
    disp('The length of the raw FID data is longer than the designated zero-filling length.')
    disp('\nOverwrite the zero-filling length with the FID length.')
    JET_STR.SR.params.ZeroFillTo = size(JET_STR.SR.data.organized_fid_data, 1);
end

% For Version 1, only do the zero-filling.
fid_data_less_processing = zeros(JET_STR.SR.params.ZeroFillTo, size(JET_STR.SR.data.organized_fid_data, 2));
fid_data_less_processing(1 : size(JET_STR.SR.data.organized_fid_data, 1), :) = JET_STR.SR.data.organized_fid_data;
JET_STR.SR.data.fid_data_less_processing = fid_data_less_processing;

% For Version 2, do the zero-filling as well as other subsequent pre-processing steps.
fid_data_more_processing = zeros(JET_STR.SR.params.ZeroFillTo, size(JET_STR.SR.data.organized_fid_data, 2));
fid_data_more_processing(1 : size(JET_STR.SR.data.organized_fid_data, 1), :) = JET_STR.SR.data.organized_fid_data;
% Spectral line broadening.
fid_data_more_processing_LB = fid_data_more_processing .* ...
    repmat((exp( - (JET_STR.SR.data.time_zeropad') * JET_STR.SR.params.LB * pi)), [1, JET_STR.SR.params.num_row]);
JET_STR.SR.data.fid_data_more_processing = fid_data_more_processing_LB;

% Time domain -> Frequency domain.
JET_STR.SR.data.spectrum_less_processing = ...
    fftshift(fft(JET_STR.SR.data.fid_data_less_processing, JET_STR.SR.params.ZeroFillTo, 1), 1);

JET_STR.SR.data.spectrum_more_processing = ...
    fftshift(fft(JET_STR.SR.data.fid_data_more_processing, JET_STR.SR.params.ZeroFillTo, 1), 1);

% Create the x-axis (frequency-axis).
JET_STR.SR.data.frequency_axis = (1 : 1 : JET_STR.SR.params.ZeroFillTo) / ...
    JET_STR.SR.params.ZeroFillTo * JET_STR.SR.params.ppm + 4.7 - JET_STR.SR.params.ppm / 2.0;

% Define the frequency bound for plotting.
JET_STR.SR.data.Display_freqbounds = ...
    JET_helper_function_find_freqbounds(JET_STR.SR.data.frequency_axis, ...
    JET_STR.SR.params.Display_upperbound, JET_STR.SR.params.Display_lowerbound);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Perform spectrum registration based on the SR procedure.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch JET_STR.SR.procedure
    case 'intraOnOff_interOnOff'
        JET_STR = JET_SR_procedure_intraOnOff_interOnOff(JET_STR);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

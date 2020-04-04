function corrected_spectra = JET_SR_apply_frequency_phase_correction_JETSR(spectra, JET_STR, correction_parameters)
%JET_SR_APPLY_FREQUENCY_PHASE_CORRECTION_JETSR Apply the calculated
%frequency and phase correction for each individual repetition.
%
% Input arguments:
% - SPECTRA: 2D spectra with the following shape:
%   (repetition) x (spectra length)
% - CORRECTION_PARAMETERS: The frequency and phase correction parameters
% for each repetition.
%
% Output arguments:
% CORRECTED_SPECTRA: frequency and phase corrected spectra of the same size
% as the input SPECTRA.
%
% Chen "Raphael" Liu (cl3760@columbia.edu) & Jia Guo (jg3400@columbia.edu),
% 02/24/2020.

if length(size(spectra)) ~= 2
    error(['Expecting 2D spectra but its shape is [', ...
        num2str(size(spectra)), ...
        ']. Cannot run JET_SR_apply_frequency_phase_correction_JETSR.'])
end

num_repetition = size(spectra, 1);

if length(size(correction_parameters)) ~= 2
    error(['Expecting 2D correction_parameters but its shape is [', ...
        num2str(size(correction_parameters)), ...
        ']. Cannot run JET_SR_apply_frequency_phase_correction_JETSR.'])
end

if size(correction_parameters, 1) ~= num_repetition
    error(['Expecting correction_parameters to have the same number', ...
        'of rows as the spectra but its shape is [', ...
        num2str(size(correction_parameters)), ...
        ']. Cannot run JET_SR_apply_frequency_phase_correction_JETSR.'])
end

% Apply the calculated frequency and phase correction.
corrected_spectra = zeros(size(spectra));

for repetition_index = 1 : num_repetition
    % Calculate the corresponding FID from the spectra and apply the
    % transform on the FID.
    spectra_reversecalculated_FID = ...
            ifft(fliplr(spectra(repetition_index, :)));
    corrected_spectra(repetition_index, :) = ...
        JET_helper_function_spectrum_deformation_complex(correction_parameters(repetition_index, :), ...
        JET_STR.SR.data.time_zeropad, 1 : length(JET_STR.SR.data.frequency_axis), ...
        spectra_reversecalculated_FID, JET_STR.SR.params.LS_mode);
end

end
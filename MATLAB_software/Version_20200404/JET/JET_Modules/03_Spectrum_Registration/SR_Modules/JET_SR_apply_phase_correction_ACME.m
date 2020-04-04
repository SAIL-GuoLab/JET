function phase_corrected_spectra = JET_SR_apply_phase_correction_ACME(spectra, phase_correction_parameters)
%JET_SR_APPLY_PHASE_CORRECTION_ACME Apply the initial coil channel-
%dependent, repetition-independent phase correction.
%
% Input arguments:
% - SPECTRA: 3D spectra with the following shape:
%   (coil channel) x (repetition) x (spectra length)
% - PHASE_CORRECTION_PARAMETERS: The frequency correction parameters
% for each coil channel.
%
% Output arguments:
% PHASE_CORRECTED_SPECTRA: phase corrected spectra of the same size as the
% input SPECTRA.
%
% Chen "Raphael" Liu (cl3760@columbia.edu) & Jia Guo (jg3400@columbia.edu),
% 02/24/2020.

if length(size(spectra)) ~= 3
    error(['Expecting 3D spectra but its shape is [', ...
        num2str(size(spectra)), ...
        ']. Cannot run JET_SR_apply_phase_correction_ACME.'])
end

num_channel = size(spectra, 1);
num_repetition = size(spectra, 2);

if length(size(phase_correction_parameters)) ~= 2
    error(['Expecting 2D phase_correction_parameters but its shape is [', ...
        num2str(size(phase_correction_parameters)), ...
        ']. Cannot run JET_SR_apply_phase_correction_ACME.'])
end

if size(phase_correction_parameters, 1) ~= num_channel
    error(['Expecting phase_correction_parameters to have the same number', ...
        'of rows as the spectra but its shape is [', ...
        num2str(size(phase_correction_parameters)), ...
        ']. Cannot run JET_SR_apply_phase_correction_ACME.'])
end

% Apply the calculated phase correction.
phase_corrected_spectra = zeros(size(spectra));

for channel_index = 1 : num_channel
    for repetition_index = 1 : num_repetition
        phase_corrected_spectra(channel_index, repetition_index, :) = ...
            dephase_fun(spectra(channel_index, repetition_index, :), ...
            - phase_correction_parameters(channel_index, :));
    end
end

end

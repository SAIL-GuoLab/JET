function frequency_corrected_spectra = JET_SR_apply_frequency_correction_ICOSHIFT(spectra, frequency_correction_parameters)
%JET_SR_APPLY_FREQUENCY_CORRECTION_ICOSHIFT Apply the initial coil channel-
%dependent, repetition-independent frequency correction.
%
% Input arguments:
% - SPECTRA: 3D spectra with the following shape:
%   (coil channel) x (repetition) x (spectra length)
% - FREQUENCY_CORRECTION_PARAMETERS: The frequency correction parameters
% for each coil channel.
%
% Output arguments:
% FREQUENCY_CORRECTED_SPECTRA: frequency corrected spectra of the same size
% as the input SPECTRA.
%
% Chen "Raphael" Liu (cl3760@columbia.edu) & Jia Guo (jg3400@columbia.edu),
% 02/24/2020.

if length(size(spectra)) ~= 3
    error(['Expecting 3D spectra but its shape is [', ...
        num2str(size(spectra)), ...
        ']. Cannot run JET_SR_apply_frequency_correction_ICOSHIFT.'])
end

num_channel = size(spectra, 1);
num_repetition = size(spectra, 2);

if length(size(frequency_correction_parameters)) ~= 2
    error(['Expecting 2D frequency_correction_parameters but its shape is [', ...
        num2str(size(frequency_correction_parameters)), ...
        ']. Cannot run JET_SR_apply_frequency_correction_ICOSHIFT.'])
end

if size(frequency_correction_parameters, 1) ~= num_channel
    error(['Expecting frequency_correction_parameters to have the same number', ...
        'of rows as the spectra but its shape is [', ...
        num2str(size(frequency_correction_parameters)), ...
        ']. Cannot run JET_SR_apply_frequency_correction_ICOSHIFT.'])
end

% Apply the calculated frequency correction.
frequency_corrected_spectra = zeros(size(spectra));

for channel_index = 1 : num_channel
    for repetition_index = 1 : num_repetition
        frequency_corrected_spectra(channel_index, repetition_index, :) = ...
            my_shift(spectra(channel_index, repetition_index, :), ...
            frequency_correction_parameters(channel_index, :));
    end
end

end
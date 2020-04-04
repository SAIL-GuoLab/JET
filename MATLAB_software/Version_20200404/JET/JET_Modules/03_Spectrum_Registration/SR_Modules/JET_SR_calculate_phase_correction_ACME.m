function phase_correction_parameters = JET_SR_calculate_phase_correction_ACME(spectra)
%JET_SR_CALCULATE_PHASE_CORRECTION_ACME Use ACME to calculate parameters
%for the initial coil channel-dependent, repetition-independent phase
%correction.
%
% Input arguments:
% - SPECTRA: 3D spectra with the following shape:
%   (coil channel) x (repetition) x (spectra length)
%
% Output arguments:
% - PHASE_CORRECTION_PARAMETERS: The frequency correction parameters
% for each coil channel.
%
% Chen "Raphael" Liu (cl3760@columbia.edu) & Jia Guo (jg3400@columbia.edu),
% 02/24/2020.

if length(size(spectra)) ~= 3
    error(['Expecting 3D spectra but its shape is [', ...
        num2str(size(spectra)), ...
        ']. Cannot run JET_SR_calculate_phase_correction_ACME.'])
end

num_channel = size(spectra, 1);

% Initialize the parameters as zeros.
phase_correction_parameters = zeros(num_channel, 1);
% For each coil channel, take the mean across all repetitions and use
% that as a high-SNR representation for that coil channel.
spectra_SepChannel_CombRep = mean(spectra, 2);
if num_channel == 1
    spectra_SepChannel_CombRep = squeeze(spectra_SepChannel_CombRep)';
else
    spectra_SepChannel_CombRep = squeeze(spectra_SepChannel_CombRep);
end

for channel_index = 1 : num_channel
    % Calculate the coil channel-dependent, repetition-independent phase
    % shift to correct zero order phase.
    [~, phase_shift] = ...
        ACME(spectra_SepChannel_CombRep(channel_index, :), 0);
    phase_shift = mod(phase_shift, 360);
    if phase_shift > 180
        phase_shift = phase_shift - 360;
    end
    phase_shift = min(phase_shift, 360 - phase_shift);
    phase_correction_parameters(channel_index, 1) = phase_shift;
end

end

function frequency_correction_parameters = JET_SR_calculate_frequency_correction_ICOSHIFT(spectra, icoshift_freqbounds)
%JET_SR_FREQUENCY_CORRECTION_ICOSHIFT Use ICOSHIFT to complete the initial
%coil channel-dependent, repetition-independent frequency correction.
%
% Input arguments:
% - SPECTRA: 3D spectra with the following shape:
%   (coil channel) x (repetition) x (spectra length)
% - ICOSHIFT_FREQBOUNDS: The frequency bounds over which we minimize the
%   differences across coil channels.
%
% Output arguments:
% - FREQUENCY_CORRECTION_PARAMETERS: The frequency correction parameters
% for each coil channel.
%
% Chen "Raphael" Liu (cl3760@columbia.edu) & Jia Guo (jg3400@columbia.edu),
% 02/24/2020.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Initial frequency correction with icoshift.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(size(spectra)) ~= 3
    error(['Expecting 3D spectra but its shape is [', ...
        num2str(size(spectra)), ...
        ']. Cannot run JET_SR_calculate_frequency_correction_ICOSHIFT.'])
end

num_channel = size(spectra, 1);

% For each coil channel, take the mean across all repetitions and use
% that as a high-SNR representation for that coil channel.
spectra_SepChannel_CombRep = mean(spectra, 2);
if num_channel == 1
    spectra_SepChannel_CombRep = squeeze(spectra_SepChannel_CombRep)';
else
    spectra_SepChannel_CombRep = squeeze(spectra_SepChannel_CombRep);
end

sqr_abs_spectra = abs(spectra_SepChannel_CombRep(:, icoshift_freqbounds)) .^ 2;

% Calculate the coil channel-dependent, repetition-independent frequency shift.
frequency_correction_parameters = zeros(num_channel, 1);
for icoshift_iter = 1 : 3
    [~, ~, icoshift_indices] = ...
        icoshift('average', sqr_abs_spectra, ...
        'whole', 'f', [2, 1, 0], ...
        1 : 1 : size(sqr_abs_spectra, 2));
    close;
    if sum(icoshift_indices) == 0
        frequency_correction_parameters(:, 1) = frequency_correction_parameters(:, 1) + icoshift_indices;
        break
    end
    frequency_correction_parameters(:, 1) = frequency_correction_parameters(:, 1) + icoshift_indices;
end

end
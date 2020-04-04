function [correction_parameters, spectra_RepSmooth_similarity, spectra_SR_template] = ...
    JET_SR_calculate_frequency_phase_correction_JETSR(spectra, JET_STR, SR_freqbounds, spectra_SR_template)
%JET_SR_FREQUENCY_PHASE_CORRECTION_JETSR Use JET SR code to complete the
%follow-up spectrum registration by correcting the frequencies and 
%zero-order phases for each individual repetition.
%We assume that by this point the spectra shall already been 
%channel-combined (if applicable).
%
% Input arguments:
% - SPECTRA: 2D spectra with the following shape:
%   (repetition) x (spectra length)
% - JET_STR: JET struct.
% - SR_FREQBOUNDS: The frequency bounds over which we minimize the
%   differences across each repetition.
%
% Output arguments:
% - CORRECTION_PARAMETERS: The correction parameters (frequency and
% zero-order phase) for each coil channel.
%
% Chen "Raphael" Liu (cl3760@columbia.edu) & Jia Guo (jg3400@columbia.edu),
% 02/24/2020.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Spectrum registration by correcting frequencies for each repetition.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(size(spectra)) ~= 2
    error(['Expecting 2D spectra but its shape is [', ...
        num2str(size(spectra)), ...
        ']. Cannot run JET_SR_calculate_frequency_phase_correction_JETSR.'])
end

if ~exist('spectra_SR_template', 'var')
    SR_template_exist = 0;
else
    SR_template_exist = 1;
end

num_repetition = size(spectra, 1);
spectra_length = size(spectra, 2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Inter-repetition smoothing.
%%%% We will use the smoothed version to calculate the spectrum
%%%% registration tranforms, and apply the transforms on the original data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gaussian_window = gausswin(num_repetition, JET_STR.SR.params.SR_RepSmooth_Alpha);
gaussian_window_normalized = gaussian_window ./ sum(gaussian_window);
gaussian_filter = @(signal) conv(signal, gaussian_window_normalized, 'same');

spectra_RepSmooth = zeros(num_repetition, spectra_length);
for frequency_index = 1 : spectra_length
    spectra_RepSmooth(:, frequency_index) = gaussian_filter(spectra(:, frequency_index));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% How do we optimize the parameters for spectrum registration?
% We define the pair [chemical_shift, zero_order_phase], and optimize the 
% resulting spectrum. That is to say, we need to find the FID from the 
% initial moving spectrum first, and then apply the 
% [chemical_shift, zero_order_phase] pair onto the FID and find the
% parameters that will yield a resulting spectrum that best resembles the
% template spectrum.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% COMPONENT 0. Initializations. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Initialize matrices to store the registration parameters
% [chemical_shift, zero_order_phase].
correction_parameters = zeros(num_repetition, 2);
% 2. Initialize the updated spectra.
spectra_updated = zeros(num_repetition, spectra_length);
% 3. Initialize similarity matrices. to represent the
% similarity between each individual spectrum and the template.
spectra_RepSmooth_similarity = zeros(num_repetition, 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for SR_iteration_index = 1 : JET_STR.SR.params.SR_iterations
    %%%%% COMPONENT 1. Prepare the spectrum templates. %%%%%%%%%%%%%%%%
    % If the SR template is already given, ignore this step. Otherwise
    % build the template using the following approach.
    if SR_template_exist == 0
        % For the first iteration, use the mean spectrum as the template.
        if SR_iteration_index == 1
            if strcmp(JET_STR.SR.params.SR_with_real_or_complex, 'real')
                spectra_SR_template = mean(real(spectra_RepSmooth), 1);
            elseif strcmp(JET_STR.SR.params.SR_with_real_or_complex, 'complex')
                spectra_SR_template = mean(spectra_RepSmooth, 1);
            end
        % For subsequent iterations, use the weighted sum of each
        % individual spectrum (with the weights proportional to the
        % similarity between each individual spectrum and the current
        % template) as the new template.
        else
            spectra_SR_weightings = ...
                repmat(spectra_RepSmooth_similarity_normalized , [1, spectra_length]);
            
            if strcmp(JET_STR.SR.params.SR_with_real_or_complex, 'real')
                spectra_SR_template = ...
                    sum(real(spectra_RepSmooth_updated) .* spectra_SR_weightings, 1);
            elseif strcmp(JET_STR.SR.params.SR_with_real_or_complex, 'complex')
                spectra_SR_template = ...
                    sum(spectra_RepSmooth_updated .* spectra_SR_weightings, 1);
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%% COMPONENT 2. Register each spectrum to the spectrum template.
    % Initialize spectral deformation parameters.
    % [chemical_shift, zero_order_phase]
    SR_params_initial = [0, 0];
    for repetition_index = 1 : num_repetition
        % Reverse-calculate the moving FID from the initial moving spectrum.
        spectra_reversecalculated_FID = ...
            ifft(fliplr(spectra(repetition_index, :)));
        spectra_RepSmooth_reversecalculated_FID = ...
            ifft(fliplr(spectra_RepSmooth(repetition_index, :)));
        
        if strcmp(JET_STR.SR.params.SR_with_real_or_complex, 'real')
            % Spectrum Registration, within On spectra.
            [singleRep_SR_params] = lsqcurvefit(@(parameters, SR_params_initial) ...
                JET_helper_function_spectrum_deformation_real(parameters, ...
                SR_params_initial, SR_freqbounds, spectra_RepSmooth_reversecalculated_FID, ...
                JET_STR.SR.params.LS_mode), SR_params_initial, JET_STR.SR.data.time_zeropad, ...
                real(spectra_SR_template(SR_freqbounds)));
        elseif strcmp(JET_STR.SR.params.SR_with_real_or_complex, 'complex')
            % Spectrum Registration, within On spectra.
            [singleRep_SR_params] = lsqcurvefit(@(parameters, SR_params_initial) ...
                JET_helper_function_spectrum_deformation_complex(parameters, ...
                SR_params_initial, SR_freqbounds, spectra_RepSmooth_reversecalculated_FID, ...
                JET_STR.SR.params.LS_mode), SR_params_initial, JET_STR.SR.data.time_zeropad, ...
                spectra_SR_template(SR_freqbounds));
        end
        
        % Update the spectrum registration parameters for the current repetition.
        correction_parameters(repetition_index, :) = ...
            singleRep_SR_params + correction_parameters(repetition_index, :);
        
        % Calculate the new spectrum after this current iteration of
        % spectrum registration. Remember to use the version without
        % inter-repetition smoothing.
        spectra_updated(repetition_index, :) = ...
            JET_helper_function_spectrum_deformation_complex(correction_parameters(repetition_index, :), ...
            JET_STR.SR.data.time_zeropad, 1 : length(JET_STR.SR.data.frequency_axis), ...
            spectra_reversecalculated_FID, JET_STR.SR.params.LS_mode);
    end
    
    % Apply inter-repetition smoothing on the newly registered spectra.
    spectra_RepSmooth_updated = zeros(num_repetition, spectra_length);
    for frequency_index = 1 : spectra_length
        spectra_RepSmooth_updated(:, frequency_index) = ...
            gaussian_filter(spectra_updated(:, frequency_index));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%% COMPONENT 3. Calculate and normalize similarity %%%%%%%%%%%%%
    % Similarity is defined between each individual spectrum and the
    % current template.
    if strcmp(JET_STR.SR.params.SR_with_real_or_complex, 'real')
        for repetition_index = 1 : num_repetition
            correlation_coefficient_matrix = ...
                corrcoef([real(spectra_SR_template(SR_freqbounds))', ...
                real(spectra_RepSmooth_updated(repetition_index, SR_freqbounds))']);
            spectra_RepSmooth_similarity(repetition_index, :) = ...
                correlation_coefficient_matrix(1, 2);
        end
    elseif strcmp(JET_STR.SR.params.SR_with_real_or_complex, 'complex')
        for repetition_index = 1 : num_repetition
            correlation_coefficient_matrix = ...
                corrcoef([spectra_SR_template(SR_freqbounds)', ...
                spectra_RepSmooth_updated(repetition_index, SR_freqbounds)']);
            spectra_RepSmooth_similarity(repetition_index, :) = ...
                abs(correlation_coefficient_matrix(1, 2));
        end
    end
    
    % Normalize the similarity matrices such that the sum is 1. This is
    % necessary for creating a new spectrum template as the weighted
    % sum of all spectra.
    spectra_RepSmooth_similarity_normalized = ...
        spectra_RepSmooth_similarity ./ sum(spectra_RepSmooth_similarity);
    
    % If this is the last iteration, we will update the templates once
    % again. Of course this does not apply if the SR template is given as
    % the input.
    if SR_template_exist == 0
        if SR_iteration_index == JET_STR.SR.params.SR_iterations
            spectra_SR_weightings = ...
                repmat(spectra_RepSmooth_similarity_normalized, [1, spectra_length]);
            
            if strcmp(JET_STR.SR.params.SR_with_real_or_complex, 'real')
                spectra_SR_template = ...
                    sum(real(spectra_RepSmooth_updated) .* spectra_SR_weightings, 1);
            elseif strcmp(JET_STR.SR.params.SR_with_real_or_complex, 'complex')
                spectra_SR_template = ...
                    sum(spectra_RepSmooth_updated .* spectra_SR_weightings, 1);
            end
        end
    end
end

end
function organized_FIDs = JET_Loader_load_Bruker_FIDs_VolumeCoil(gabafile, JET_STR)
% JET_LOADER_LOAD_BRUKER_FIDS_VOLUMECOIL A helper function that reads the
% FIDs, specifically for vendor == Bruker and coil_type == volume_coil.
% 
% Chen "Raphael" Liu (cl3760@columbia.edu) and Jia Guo (jg3400@columbia.edu), 02/24/2020.

disp('>>> Load the FIDs from Bruker...')

file_ID = fopen(gabafile, 'r');

if (file_ID < 0)
    errordlg('File not found.', 'File Error');
    beep;
else
    file_ID = fopen(gabafile, 'r');
    indata = fread(file_ID, inf, 'long');
    fclose(file_ID);

    size_FID = length(indata) / 2;
    indata_reshaped = reshape(indata, 2, size_FID);
    complex_data = complex(indata_reshaped(1, :), indata_reshaped(2, :));
    FID_flat = reshape(complex_data, size_FID, 1);

    % Initialize the FID data matrix.
    raw_FIDs = zeros(JET_STR.SR.params.SpecLength - JET_STR.SR.params.num_compensatepoints, size(FID_flat, 1) ./ JET_STR.SR.params.SpecLength);
    for repetition_index = 1 : size(FID_flat, 1) / JET_STR.SR.params.SpecLength
        FID_current_repetition = FID_flat(1 + (repetition_index - 1) * JET_STR.SR.params.SpecLength : ...
            repetition_index * JET_STR.SR.params.SpecLength);
        
        % Throw away the first several noisy FID points per repetition.
        if (JET_STR.SR.params.num_compensatepoints > 0)
            FID_current_repetition = FID_current_repetition(JET_STR.SR.params.num_compensatepoints + 1 : end);
        end
        
        % Append the current repetition FID to the entire FID matrix.
        raw_FIDs(:, repetition_index) = FID_current_repetition';
    end

    switch JET_STR.SR.params.ONOFForder
        case 'onfirst'
            ON_OFF_sequence = ~rem(idivide(int16(1 : size(raw_FIDs, 2)) - 1, JET_STR.SR.params.num_average), 2);
        case 'offfirst'
            ON_OFF_sequence = rem(idivide(int16(1 : size(raw_FIDs, 2)) - 1, JET_STR.SR.params.num_average), 2);
    end

    FIDs_On = raw_FIDs(:, (ON_OFF_sequence == 1));
    FIDs_Off = raw_FIDs(:, (ON_OFF_sequence == 0));
    
    FIDs_On_Off_Interleave = zeros(size(FIDs_On, 1), 2 * size(FIDs_On, 2));
    FIDs_On_Off_Interleave(:, 1:2:end) = FIDs_On;
    FIDs_On_Off_Interleave(:, 2:2:end) = FIDs_Off;
    
    organized_FIDs = FIDs_On_Off_Interleave;
    
    disp('...FID Loaded! Re-organized as On-Off interleave spectra.')
end
end
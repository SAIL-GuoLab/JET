function organized_FIDs = JET_Loader_load_Bruker_FIDs_MultiArraySurfaceCoil(gabafile, JET_STR)
% JET_LOADER_LOAD_BRUKER_FIDS_VOLUMECOIL A helper function that reads the
% FIDs, specifically for vendor == Bruker and coil_type ==
% multi_array_surface_coil.
%
% Chen "Raphael" Liu (cl3760@columbia.edu) and Jia Guo (jg3400@columbia.edu), 02/24/2020.

disp('>>> Load the FIDs from Bruker Multi-Array Surface Coil...')

file_ID = fopen(gabafile, 'r');

if (file_ID < 0)
    errordlg('File not found.', 'File Error');
    beep;

else
    % Read the FIDs from the file.
    file_ID = fopen(gabafile, 'r');
    indata = fread(file_ID, inf, 'long');
    fclose(file_ID);

    % Extract the FIDs and seperate them into correct number of rows and
    % columns ((points per spectrum) x (number of spectra)). This is stored
    % in the variable "raw_FIDs".
    size_FID = length(indata) / 2;
    indata_reshaped = reshape(indata, 2, size_FID);
    complex_data = complex(indata_reshaped(1, :), indata_reshaped(2, :));
    complex_data = reshape(complex_data, size_FID, 1);
    reshaped_data = reshape(complex_data, JET_STR.SR.params.SpecLength, size_FID / JET_STR.SR.params.SpecLength);
    raw_FIDs = reshaped_data(JET_STR.SR.params.num_compensatepoints + 1 : end, :);

    % Reorganize the FIDs into the following shape:
    % (coil channel x repetition) x (FID length)
    %switch JET_STR.SR.params.ONOFForder
    %    case 'onfirst'
    %        ON_OFF_sequence = ~rem(idivide(int16(1 : size(raw_FIDs, 2)) - 1, 1), 2);
    %    case 'offfirst'
    %        ON_OFF_sequence = rem(idivide(int16(1 : size(raw_FIDs, 2)) - 1, 1), 2);
    %end
    
    %FIDs_On = raw_FIDs(:, (ON_OFF_sequence == 1));
    %FIDs_Off = raw_FIDs(:, (ON_OFF_sequence == 0));
    
    %FIDs_On_Off_Interleave = zeros(size(FIDs_On, 1), 2 * size(FIDs_On, 2));
    %FIDs_On_Off_Interleave(:, 1:2:end) = FIDs_On;
    %FIDs_On_Off_Interleave(:, 2:2:end) = FIDs_Off;
    
    %organized_FIDs = FIDs_On_Off_Interleave;
    organized_FIDs = raw_FIDs;
    
    disp('...FID Loaded! Re-organized as On-Off interleave spectra.')
end
end
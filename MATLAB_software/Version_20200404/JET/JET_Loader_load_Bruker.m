function JET_STR = JET_Loader_load_Bruker(JET_STR, current_subject_folder_index)
%JET_LOADER_LOAD_BRUKER A helper function to facilitate 
%JET_Loader_update_params, which itself is part of the first step of the 
%JET pipeline. This function loads the data when the vendor is detected or
%defined as 'Bruker'.
%
% The objectives of this functions include:
% 1) read the FIDs as well as other essential parameters
% 2) correctly organize the FIDs into On-Off interleave FIDs
% 3) update essential parameters to JET_STR
%
% Input arguments
% - JET_STR : The JET struct for the study.
% - current_subject_folder_index : The index of the current folder for the
% current subject. 1 = first folder, 2 = second folder, etc.
%
% Output arguments
% - JET_STR : The modified JET struct for the study.
%
% Chen "Raphael" Liu (cl3760@columbia.edu) and Jia Guo (jg3400@columbia.edu), 02/24/2020.

% First we want to confirm that the vendor is indeed Bruker.
assert(strcmp(JET_STR.Loader.vendor, 'Bruker'), ...
    'Well the function "JET_Loader_update_params" seems to have a bug.')

% Display the progress.
disp('>>> Loading data from the Bruker vendor...')

% Find the foldername of the current subject, and search for subfolders
% that hold the scans of this subject.
current_subject_foldername = JET_STR.Loader.subject_foldername_all(current_subject_folder_index).name;
current_subject_path = strcat(JET_STR.Loader.study_pathname, current_subject_foldername);
current_scan_foldername_all = dir(current_subject_path);
% We need to make sure that we remove '.' (the current directory) and
% '..' (the parent directory) from the search results.
current_scan_foldername_all = current_scan_foldername_all(...
    ~ismember({current_scan_foldername_all.name}, {'.', '..'}));

% Iterate over the scans of the current subject.
for scan_index = 1 : length(current_scan_foldername_all)
    current_scan_foldername = current_scan_foldername_all(scan_index).name;
    current_scan_path = strcat(JET_STR.Loader.study_pathname, current_subject_foldername, '/', current_scan_foldername);

    % Update the JET struct.
    JET_STR.Loader.MEGAPRESS_filelist{scan_index} = strcat(current_scan_path, '/rawdata.job0');
    %JET_STR.Loader.waterfile_filelist{scan_index} = strcat(current_scan_path, '/rawdata.job1');
    JET_STR.Loader.methodnamelist{scan_index} = strcat(current_scan_path, '/method');
end

for scan_index = 1 : length(JET_STR.Loader.MEGAPRESS_filelist)
    
    gabafile = JET_STR.Loader.MEGAPRESS_filelist{scan_index};
    methodname = JET_STR.Loader.methodnamelist{scan_index};

    % Try to read the header file and update the initial parameters.
    try
        file_ID = fopen(methodname);
        sparheader_struct = textscan(file_ID, '%s');
        sparheader = sparheader_struct{1};

        TEsparidx = find(cellfun('length', regexp(sparheader, 'PVM_EchoTime')) == 1); %%% Bruker Default
        TEsparidx = TEsparidx(1);
        TE_strcell = sparheader{TEsparidx};
        TE_strc = char(TE_strcell);
        JET_STR.SR.params.TE = str2double(TE_strc(17:18));

        TRsparidx = find(cellfun('length', regexp(sparheader, 'PVM_RepetitionTime')) == 1); %%% Bruker Default
        TRsparidx = TRsparidx(1);
        TR_strcell = sparheader{TRsparidx};
        TR_strc = char(TR_strcell);
        JET_STR.SR.params.TR = str2double(TR_strc(23:26));

        NumRepetitionsparidx = find(cellfun('length', regexp(sparheader, 'NRepetitions')) == 1); %%% Bruker Default
        NumRepetitionsparidx = NumRepetitionsparidx(1);
        NumRepetition_strcell = sparheader{NumRepetitionsparidx};
        NumRepetition_strc = char(NumRepetition_strcell);
        NumRepetition = NumRepetition_strc(21:end);

        NumAveragesparidx = find(cellfun('length',regexp(sparheader,'NAverage')) == 1); %%% Bruker Default
        NumAveragesparidx = NumAveragesparidx(1);
        NumAverage_strcell = sparheader{NumAveragesparidx};
        NumAverage_str = char(NumAverage_strcell);
        NumAverage = NumAverage_str(18:end);
   
        SWsparidx = find(cellfun('length', regexp(sparheader, 'PVM_DigSw')) == 1); %%% Bruker Default
        SWsparidx = SWsparidx(1);
        SW_strcell = sparheader{SWsparidx};
        SW_strc = char(SW_strcell);
        SW = SW_strc(14:end);
        JET_STR.SR.params.sw = str2double(SW);

        PPMsparidx = find(cellfun('length', regexp(sparheader, 'PVM_SpecSW=')) == 1); %%% Bruker Default
        PPMsparidx = PPMsparidx(1) + 3;
        PPM_strcell = sparheader{PPMsparidx};
        PPM_strc = char(PPM_strcell);
        PPM = PPM_strc(1:end);
        JET_STR.SR.params.ppm = str2double(PPM);

        NPsparidx = find(cellfun('length', regexp(sparheader, 'PVM_DigNp')) == 1); %%% Bruker Default
        NPsparidx = NPsparidx(1);
        NP_strcell = sparheader{NPsparidx};
        NP_strc = char(NP_strcell);
        NP = NP_strc(14:end);
        JET_STR.SR.params.SpecLength = str2double(NP);
        JET_STR.SR.params.num_point = str2double(NP);

        AcqTimesparidx = find(cellfun('length', regexp(sparheader, 'SpecAcquisitionTime')) == 1); %%% Bruker Default
        AcqTimesparidx = AcqTimesparidx(1);
        AcqTime_strcell = sparheader{AcqTimesparidx};
        AcqTime_strc = char(AcqTime_strcell);
        AcqTime = AcqTime_strc(28:end);
        JET_STR.SR.params.AcqTime = str2double(AcqTime);

        DigShiftsparidx = find(cellfun('length', regexp(sparheader, 'PVM_DigShift')) == 1); %%% Bruker Default
        DigShiftsparidx = DigShiftsparidx(1);
        DigShift_strcell = sparheader{DigShiftsparidx};
        DigShift_strc = char(DigShift_strcell);
        DigShift = DigShift_strc(17:end);
        JET_STR.SR.params.num_compensatepoints = str2double(DigShift);

        VoxArrSizeidx = find(cellfun('length', regexp(sparheader, 'PVM_VoxArrSize')) == 1); %%% Bruker Default
        VoxArrSizeidx = VoxArrSizeidx(1);
        JET_STR.SR.params.voxsize(1,2) = str2double(char(sparheader{VoxArrSizeidx+4})); %'ap_size'???
        JET_STR.SR.params.voxsize(1,1) = str2double(char(sparheader{VoxArrSizeidx+5})); %'lr_size'???
        JET_STR.SR.params.voxsize(1,3) = str2double(char(sparheader{VoxArrSizeidx+6})); %'cc_size'???
        JET_STR.SR.params.num_point = JET_STR.SR.params.num_point - JET_STR.SR.params.num_compensatepoints;
        
        JET_STR.SR.params.num_repetition = str2double(NumRepetition);
        JET_STR.SR.params.num_average = str2double(NumAverage);
        
    catch
    end

    % Load all the FIDs, length x averages.
    if scan_index == 1
        JET_STR.Loader.data.organized_fid_data = JET_Loader_load_Bruker_FIDs(gabafile, JET_STR);
    else
        JET_STR.Loader.data.organized_fid_data = [JET_STR.Loader.data.organized_fid_data, JET_Loader_load_Bruker_FIDs(gabafile, JET_STR)];
    end
    
    JET_STR.SR.params.num_row = size(JET_STR.Loader.data.organized_fid_data, 2);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Double check if the FIDs extracted look reasonable.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Spectrum = fftshift(fft(JET_STR.Loader.data.organized_fid_data, [], 1), 1);
Spectrum_On = Spectrum(:, 1:2:end);
Spectrum_Off = Spectrum(:, 2:2:end);
frequency_axis = (1 : 1 : size(Spectrum, 1)) / ...
    size(Spectrum, 1) * JET_STR.SR.params.ppm ...
    + JET_STR.SR.params.water_ppm - JET_STR.SR.params.ppm / 2.0;

figure(1);
imagesc(real([Spectrum_On'; Spectrum_Off']));
colorbar();
title('Raw On (upper half) and Off (lower half) Spectra');

figure(2);
subplot(2, 1, 1)
plot(frequency_axis, real(Spectrum)); title('Raw Real Spectra');
set(gca, 'Xdir', 'reverse')

subplot(2, 1, 2)
plot(frequency_axis, mean(real(Spectrum_On), 2));
hold on;
plot(frequency_axis, mean(real(Spectrum_Off), 2));
title('Raw On and Off Spectra Mean');
legend('On', 'Off');
set(gca, 'Xdir', 'reverse')

drawnow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('...Data loaded from Bruker!')

end
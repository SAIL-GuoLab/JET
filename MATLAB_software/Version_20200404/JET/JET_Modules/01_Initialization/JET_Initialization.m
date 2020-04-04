function JET_STR = JET_Initialization(JET_STR)
%JET_LOADER Preparation step (even before the first step) in the JET
% pipeline: initialization.
%
% Input arguments
% - JET_STR : The JET struct for the study.
%
% Output arguments
% - JET_STR : The modified JET struct for the study.
%
% Chen "Raphael" Liu (cl3760@columbia.edu) & Jia Guo (jg3400@columbia.edu),
% 02/24/2020.

% Display the progress.
disp('>>>>>> JET Initialization begin...');
% Initialize output metabolites.
JET_STR.Report.metabolites.GABA_Amplitude= zeros(JET_STR.Report.metabolites.num_samples, 1);
JET_STR.Report.metabolites.Glx_Amplitude = zeros(JET_STR.Report.metabolites.num_samples, 1);
JET_STR.Report.metabolites.Cr_Amplitude = zeros(JET_STR.Report.metabolites.num_samples, 1);
JET_STR.Report.metabolites.Cho_Amplitude = zeros(JET_STR.Report.metabolites.num_samples, 1);
JET_STR.Report.metabolites.NAA_Amplitude = zeros(JET_STR.Report.metabolites.num_samples, 1);
JET_STR.Report.metabolites.Glx_GOF = zeros(JET_STR.Report.metabolites.num_samples, 1);
JET_STR.Report.metabolites.GABA_GOF = zeros(JET_STR.Report.metabolites.num_samples, 1);
JET_STR.Report.metabolites.Cr_GOF = zeros(JET_STR.Report.metabolites.num_samples, 1);
JET_STR.Report.metabolites.Cho_GOF = zeros(JET_STR.Report.metabolites.num_samples, 1);
JET_STR.Report.metabolites.NAA_GOF = zeros(JET_STR.Report.metabolites.num_samples, 1);
JET_STR.Report.metabolites.Cr_FWHM = zeros(JET_STR.Report.metabolites.num_samples, 1);
JET_STR.Report.metabolites.NSDR_GABA = zeros(JET_STR.Report.metabolites.num_samples, 1);
JET_STR.Report.metabolites.NSDR_Glx = zeros(JET_STR.Report.metabolites.num_samples, 1);

% Search for all the subject folders of the study.
% We need to make sure that we remove '.' (the current directory) and 
% '..' (the parent directory) from the search results.
JET_STR.Loader.subject_foldername_all = dir(JET_STR.Loader.study_pathname);
JET_STR.Loader.subject_foldername_all = ...
    JET_STR.Loader.subject_foldername_all(...
    ~ismember({JET_STR.Loader.subject_foldername_all.name}, {'.', '..'}));

% Build a table to store the future report.
% Create a directory for the table unless the directory already exists.
if (exist(JET_STR.Report.table.foldername, 'dir') ~= 7)
    eval(['mkdir ', JET_STR.Report.table.foldername]);
end

% Define the header of the table.
JET_STR.Report.table.header = ...
    ["Subject_Name", "GABA_Amplitude", "Glx_Amplitude", "Cr_Amplitude", "Cho_Amplitude", ...
    "NAA_Amplitude", "GABA_GOF", "Glx_GOF", "Cr_GOF", "Cho_GOF", "NAA_GOF", "Cr_FWHM", ...
    "NSDR_GABA", "NSDR_Glx"];

% Create a file placeholder for the table.
JET_STR.Report.table.file_ID = fopen(JET_STR.Report.table.filename, 'w');
% special note: [~] is to prevent the system from assigning a value to the
% variable 'ans', which is a complete nonsense.
[~] = fprintf(JET_STR.Report.table.file_ID, '%s,', JET_STR.Report.table.header);
[~] = fclose(JET_STR.Report.table.file_ID);

% Delete the file_ID, which is useless beyond this point.
JET_STR.Report.table = rmfield(JET_STR.Report.table, 'file_ID');

disp('... JET initialization done!');
end
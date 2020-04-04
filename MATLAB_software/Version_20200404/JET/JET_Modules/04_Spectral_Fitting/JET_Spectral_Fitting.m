function JET_STR = JET_Spectral_Fitting(JET_STR, current_subject_folder_index)
%JET_Spectral_Fitting Third step in the JET pipeline: Perform 
% spectral fitting.
%
% Input arguments
% - JET_STR : The JET struct for the study.
% - current_subject_folder_index : The index of the current folder for the
% current subject. 1 = first folder, 2 = second folder, etc.
%
% Output arguments
% - JET_STR : The modified JET struct for the study.
%
% Chen "Raphael" Liu (cl3760@columbia.edu) & Jia Guo (jg3400@columbia.edu),
% 02/24/2020.

%% set up basic parameters
R2 = JET_STR.SF.params.LB_metabolites; % Hz
ZF = JET_STR.SF.params.ZeroFillTo;
time_zeropad = (1 : 1 : ZF) / JET_STR.SR.params.sw;

%%  STEP A - load simulated basis set: GABA, Glu, Gln, Cr, Cho, NAA ...
% load GABA basis set (Off)
GABA_basis_file = load('./NMRWizard_SimulatedBasisSet/Off/GABAH4');
GABA_basis_fid = GABA_basis_file(:, 1) + 1i .* GABA_basis_file(:, 2);
GABA_basis_fid_off = GABA_basis_fid .* exp(-(time_zeropad') * R2 * pi);

% load Glu basis set (Off)
Glu_basis_file = load('./NMRWizard_SimulatedBasisSet/Off/GlutamateH2');
Glu_basis_fid = Glu_basis_file(:, 1) + 1i .* Glu_basis_file(:, 2);
Glu_basis_fid_off = Glu_basis_fid .* exp(-(time_zeropad') * R2 * pi);

% load Gln basis set (Off)
Gln_basis_file = load('./NMRWizard_SimulatedBasisSet/Off/GlutamineH2');
Gln_basis_fid = Gln_basis_file(:, 1) + 1i .* Gln_basis_file(:, 2);
Gln_basis_fid_off = Gln_basis_fid .* exp(-(time_zeropad') * R2 * pi);

% load 2HG basis set (Off)
HG_basis_file = load('./NMRWizard_SimulatedBasisSet/Off/2HGH2');
HG_basis_fid = HG_basis_file(:, 1) + 1i .* HG_basis_file(:, 2);
HG_basis_fid_off = HG_basis_fid .* exp(-(time_zeropad') * R2 * pi);

% load GABA basis set (On)
GABA_basis_file = load('./NMRWizard_SimulatedBasisSet/On/GABAH4');
GABA_basis_fid = GABA_basis_file(:, 1) + 1i .* GABA_basis_file(:, 2);
GABA_basis_fid_on = GABA_basis_fid .* exp(-(time_zeropad') * R2 * pi);

% load Glu basis set (On)
Glu_basis_file = load('./NMRWizard_SimulatedBasisSet/On/GlutamateH2');
Glu_basis_fid = Glu_basis_file(:, 1) + 1i .* Glu_basis_file(:, 2);
Glu_basis_fid_on = Glu_basis_fid .* exp(-(time_zeropad') * R2 * pi);

% load Gln basis set (On)
Gln_basis_file = load('./NMRWizard_SimulatedBasisSet/On/GlutamineH2');
Gln_basis_fid = Gln_basis_file(:, 1) + 1i .* Gln_basis_file(:, 2);
Gln_basis_fid_on = Gln_basis_fid .* exp(-(time_zeropad') * R2 * pi);

% load 2HG basis set (On)
%HG_basis_file = load('./NMRWizard_SimulatedBasisSet/On/2HGH2');
%HG_basis_fid = HG_basis_file(:, 1) + 1i .* HG_basis_file(:, 2);
%HG_basis_fid_on = HG_basis_fid .* exp(-(time') * R2 * pi);

% create the Diff basis set for Glx, Gaba, 2HG ...
% merge Glu and Gln to create the Glx, Glu:Gln=8:2
Glx_basis_fid_on = 0.8 * Glu_basis_fid_on + 0.2 .* Gln_basis_fid_on;
Glx_basis_fid_off = 0.8 * Glu_basis_fid_off + 0.2 .* Gln_basis_fid_off;

% load Cr basis set (Off)
Cr_basis_file = load('./NMRWizard_SimulatedBasisSet/Off/Creatine');%H4
Cr_basis_fid = Cr_basis_file(:, 1) + 1i .* Cr_basis_file(:, 2);
Cr_basis_fid = Cr_basis_fid .* exp(-(time_zeropad') * R2 * pi);

% load Cho basis set (Off)
Cho_basis_file = load('./NMRWizard_SimulatedBasisSet/Off/Choline');%H3
Cho_basis_fid = Cho_basis_file(:, 1) + 1i .* Cho_basis_file(:, 2); % fid, complex number
Cho_basis_fid = Cho_basis_fid .* exp(-(time_zeropad') * R2 * pi);

% load NAA basis set (Off)
NAA_basis_file = load('./NMRWizard_SimulatedBasisSet/Off/NAA');
NAA_basis_fid = NAA_basis_file(:, 1) + 1i .* NAA_basis_file(:, 2);
NAA_basis_fid = NAA_basis_fid .* exp(-(time_zeropad') * R2 * pi);

%% normalize the basis set for fitting:  GABA, Glx, 2HG, Cr, Cho, NAA ...
% GABA, Glx, 2HG have been normalized already
Glx_basis_fid_diff = Glx_basis_fid_off - Glx_basis_fid_on;
GABA_basis_fid_diff = GABA_basis_fid_off - GABA_basis_fid_on;
%HG_basis_fid_diff = HG_basis_fid_off - HG_basis_fid_on;

Glx_basis_fid_diff = Glx_basis_fid_diff ./ max(abs(Glx_basis_fid_diff));
GABA_basis_fid_diff = GABA_basis_fid_diff ./ max(abs(GABA_basis_fid_diff));
%HG_basis_fid_diff = HG_basis_fid_diff ./ max(abs(HG_basis_fid_diff));

% reference metabolite
Cr_basis_fid = Cr_basis_fid ./ max(abs(Cr_basis_fid));
Cho_basis_fid = Cho_basis_fid ./ max(abs(Cho_basis_fid));
NAA_basis_fid = NAA_basis_fid ./ max(abs(NAA_basis_fid));

% calculate NAA basis set (Diff)
%NAA_basis_fid_diff = NAA_basis_fid .* (-1);

%% %%%%%%%%%%%%%%%%%%%%%%%% STEP B - Fitting of Diff Spectra %%%%%%%%%%%%%%%%%%%%%
% set up some parameters
DiffSpecData = fftshift(fft(ifft(fftshift(JET_STR.SR.data.diff_average_SR)) .* ...
    (exp( - (JET_STR.SR.data.time_zeropad) * JET_STR.SF.params.LB * pi))));

%% A. SQSBS - for Glx and GABA fitting
Glx_freqbounds = ...
JET_helper_function_find_freqbounds(JET_STR.SR.data.frequency_axis,...
JET_STR.SF.params.Glx_upperbound, JET_STR.SF.params.Glx_lowerbound);

GABA_freqbounds = ...
JET_helper_function_find_freqbounds(JET_STR.SR.data.frequency_axis,...
JET_STR.SF.params.GABA_upperbound, JET_STR.SF.params.GABA_lowerbound);

extra_freqbounds1 = ...
JET_helper_function_find_freqbounds(JET_STR.SR.data.frequency_axis,...
JET_STR.SF.params.extra_upperbound1, JET_STR.SF.params.extra_lowerbound1);

extra_freqbounds2 = ...
JET_helper_function_find_freqbounds(JET_STR.SR.data.frequency_axis,...
JET_STR.SF.params.extra_upperbound2, JET_STR.SF.params.extra_lowerbound2);

Glx_fitting_freqbounds = [extra_freqbounds1, Glx_freqbounds, extra_freqbounds2];
GABA_fitting_freqbounds = [extra_freqbounds1, GABA_freqbounds, extra_freqbounds2];

DiffShow_freqbounds = ...
JET_helper_function_find_freqbounds(JET_STR.SR.data.frequency_axis,...
JET_STR.SF.params.Display_upperbound, JET_STR.SF.params.Display_lowerbound);

%% Glx fit, initilization
%Glx_target_spec = DiffSpecData(Glx_fitting_freqbounds);

FID_basis = [];
FID_basis(1, :) = Glx_basis_fid_diff;
% FID_basis(3,:)=NAA_basis_fid_diff;

% check the initilization
Freq = Glx_fitting_freqbounds;

%initx = [25 30 0.14 0;32 33 0.13 0;150 70 0.16 0;0 0 0 0]; % [a_GABA R2 chemicalshift_GABA]
LS_mode = 3;

if LS_mode == 3
    initx = [1, 10, 0, 0.2, 0; 0, 0, 0, 0, 0]; % [a_GABA R2 chemicalshift_GABA]
else
    initx = [1, 10, 0.2, 0; 0, 0, 0, 0]; % [a_GABA R2 chemicalshift_GABA]
end

FitSpec = SQSBS(initx, time_zeropad, Freq, FID_basis, LS_mode);
figure();   
subplot(2, 1, 1);
plot(JET_STR.SR.data.frequency_axis(Freq), real(FitSpec));
hold on; 
plot(JET_STR.SR.data.frequency_axis(Freq), real(DiffSpecData(Freq)));
set(gca, 'XDir', 'reverse');
title('Glx Initialize');
hold off;

%% fitting
%options = optimset('lsqcurvefit');
%options = optimset(options, 'Display', 'off', 'TolFun', 1e-15, 'Tolx', 1e-10, 'MaxIter', 1e10);

%******************************************************************************************
% the core code for the least square fitting --- needs further development and optimization
[FitParams] = lsqcurvefit(@(parameters, inputx) ...
    SQSBS(parameters, inputx, Freq, FID_basis, LS_mode), initx, time_zeropad, ...
    real(DiffSpecData(Freq)));
%******************************************************************************************

% fitting result plot
Freq = Glx_freqbounds;
FitSpec = SQSBS(FitParams, time_zeropad, Freq, FID_basis, LS_mode);
Residual_FitGlx = real(DiffSpecData(Freq))' - real(FitSpec)';
JET_STR.SF.fit.FitSpec_Diff_Glx = FitSpec;
JET_STR.SF.fit.Glx_GOF = sum(Residual_FitGlx(1 : end) .^ 2) / ...
    sum(real(DiffSpecData(Glx_freqbounds(1) : Glx_freqbounds(end))) .^2 );
JET_STR.SF.fit.NSDR_Glx = std(Residual_FitGlx(1 : end));

figure()
subplot(2, 2, 1:2)
plot(JET_STR.SR.data.frequency_axis(DiffShow_freqbounds), DiffSpecData(DiffShow_freqbounds));
hold on;
plot(JET_STR.SR.data.frequency_axis(Freq), ...
    [real(FitSpec)', real(DiffSpecData(Freq))' - real(FitSpec)']);
%plot(JET_STR.SR.data.frequency_axis(Freq),real(FitSpec)');
hold on;
axis([JET_STR.SF.params.Display_lowerbound, JET_STR.SF.params.Display_upperbound, ...
    1.1 * min(real(DiffSpecData(DiffShow_freqbounds))), ...
    1.2 * max(real(DiffSpecData(DiffShow_freqbounds)))]);
set(gca, 'XDir', 'reverse');

Glx_Amplitude = FitParams(1, 1); 
%Glx_R2a = FitParams(1, 2);
%Glx_R2b = FitParams(1, 3);
%Glx_ChemShift = FitParams(1, 4);
%if size(FitParams, 2) == 5
%    Glx_Phase = FitParams(1, 5);
%elseif size(FitParams, 2) == 4
%    Glx_Phase = 0;
%end

%% GABA fit, initilization
%GABA_target_spec = DiffSpecData(GABA_fitting_freqbounds);

FID_basis = [];
FID_basis(1, :) = GABA_basis_fid_diff;
% FID_basis(3,:)=NAA_basis_fid_diff;

% check the initilization
if LS_mode == 3
    initx = [1, 10, 0, 0.2, 0; 0, 0, 0, 0, 0]; % [a_GABA R2 chemicalshift_GABA]
else
    initx = [1, 10, 0.2, 0; 0, 0, 0, 0]; % [a_GABA R2 chemicalshift_GABA]
end

Freq = GABA_fitting_freqbounds;
FitSpec = SQSBS(initx, time_zeropad, Freq, FID_basis, LS_mode);

figure();
subplot(2, 1, 1);
plot(JET_STR.SR.data.frequency_axis(Freq), real(FitSpec));
hold on; 
plot(JET_STR.SR.data.frequency_axis(Freq), real(DiffSpecData(Freq)));
set(gca, 'XDir', 'reverse');
title('GABA Initialize');
hold off;

%% GABA fitting
%options = optimset('lsqcurvefit');
%options = optimset(options,'Display','off','TolFun',1e-15,'Tolx',1e-10,'MaxIter',1e10);

%******************************************************************************************
% the core code for the least square fitting --- needs further development and optimization
[FitParams] = lsqcurvefit(@(parameters, inputx) ...
    SQSBS(parameters, inputx, Freq, FID_basis, LS_mode), initx, time_zeropad, ...
    real(DiffSpecData(Freq)));
%******************************************************************************************

% fitting result plot
Freq = GABA_freqbounds;
FitSpec = SQSBS(FitParams, time_zeropad, Freq, FID_basis, LS_mode);
Residual_FitGABA = real(DiffSpecData(Freq))' - real(FitSpec)';
JET_STR.SF.fit.FitSpec_Diff_GABA = FitSpec;

JET_STR.SF.fit.GABA_GOF = sum(Residual_FitGABA(1 : end).^2) / ...
    sum(real(DiffSpecData(GABA_freqbounds(1) : GABA_freqbounds(end))) .^2 );
JET_STR.SF.fit.NSDR_GABA = std(Residual_FitGABA(1 : end));

figure()
subplot(2, 2, 1:2)
hold on;
plot(JET_STR.SR.data.frequency_axis(Freq), ...
    [real(FitSpec)', real(DiffSpecData(Freq))' - real(FitSpec)']);
%plot(JET_STR.SR.data.frequency_axis(Freq),real(FitSpec)');
%legend('Data', 'Glx Fit', 'Glx Residual', 'GABA Fit', 'GABA Residual');
axis([JET_STR.SF.params.Display_lowerbound, JET_STR.SF.params.Display_upperbound, ...
    1.1 * min(real(DiffSpecData(DiffShow_freqbounds))), ...
    1.2 * max(real(DiffSpecData(DiffShow_freqbounds)))]);
set(gca, 'XDir', 'reverse');
title('Glx and GABA Fitted');
hold off;

GABA_Amplitude = FitParams(1, 1); 

if size(FitParams, 2) == 5
    GABA_R2a = FitParams(1, 2);
    %GABA_R2b = FitParams(1, 3);
    %GABA_ChemShift = FitParams(1, 4);
    %GABA_Phase = FitParams(1, 5);
elseif size(FitParams, 2) == 4
    GABA_R2a = FitParams(1, 2);
    %GABA_ChemShift = FitParams(1, 3);
    %GABA_Phase=FitParams(1, 4);
elseif size(FitParams, 2) == 3
    GABA_R2a = FitParams(1, 2);
    %GABA_ChemShift = FitParams(1, 3);
    %GABA_Phase = 0;
end

%% %%%%%%%%%%%%%%%%%%%%%%% Fitting of Off Spectra == LCModel of PRESS Spectra %%%%%%%%%%%%%%%%%%%%%
% prepare some parameters
OffSpecData = fftshift(fft(ifft(fftshift(JET_STR.SR.data.off_average_SR)) .* ...
    (exp( - (JET_STR.SR.data.time_zeropad) * JET_STR.SF.params.LB * pi))));

CrChoNAA_freqbounds = ...
JET_helper_function_find_freqbounds(JET_STR.SR.data.frequency_axis,...
JET_STR.SF.params.CrChoNAA_upperbound, JET_STR.SF.params.CrChoNAA_lowerbound);

%% SQSBS - for Cr&Cho&NAA fitting
%metabolite_number = 3;

Cr_lowerbound = JET_STR.SF.params.Cr_lowerbound;
Cr_upperbound = JET_STR.SF.params.Cr_upperbound;

Cho_lowerbound = JET_STR.SF.params.Cho_lowerbound;
Cho_upperbound = JET_STR.SF.params.Cho_upperbound;

NAA_lowerbound = JET_STR.SF.params.NAA_lowerbound;
NAA_upperbound = JET_STR.SF.params.NAA_upperbound;

FID_basis = [];
FID_basis(1, :) = Cr_basis_fid;
FID_basis(2, :) = Cho_basis_fid;
FID_basis(3, :) = NAA_basis_fid;

%% check the initialization
Freq = [extra_freqbounds2, CrChoNAA_freqbounds, extra_freqbounds2];

initx = [GABA_Amplitude * 5, GABA_R2a, 0, 0.25, 0; ...
    GABA_Amplitude * 2, GABA_R2a, 0, 0.25, 0; ...
    GABA_Amplitude * 5, GABA_R2a, 0, 0.35, 0;...
    0, 0, 0, 0, 0];

LS_mode = 3; % 1-Lorenzian, 2-gaussian
FitSpec = SQSBS(initx, time_zeropad, Freq, FID_basis, LS_mode);

figure();
subplot(2,1,2);
plot(JET_STR.SR.data.frequency_axis(Freq), real(FitSpec));
hold on;
plot(JET_STR.SR.data.frequency_axis(Freq), real(OffSpecData(Freq)));
set(gca, 'XDir', 'reverse');
title('Cho Cr NAA Initialize');

%% fitting the model
%options = optimset('lsqcurvefit');
%options = optimset(options,'Display','off','TolFun',1e-15,'Tolx',1e-10,'MaxIter',1e10);

%******************************************************************************************
% the core code for the least square fitting --- needs further development and optimization
[FitParams, resnorm, resid, exitflag] = ...
    lsqcurvefit(@(parameters, inputx) ...
    SQSBS(parameters, inputx, Freq, FID_basis, LS_mode), initx, time_zeropad, ...
    real(OffSpecData(Freq)));
%******************************************************************************************

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Cr_FWHM of each scan
JET_STR.SF.fit.Cr_FWHM_temp = FitParams(1, 2) / (2 * pi) + 5 + R2 - JET_STR.SR.params.LB;

%% fitting result plot
Freq = CrChoNAA_freqbounds;
FitSpec = SQSBS(FitParams, time_zeropad, Freq, FID_basis, LS_mode);
JET_STR.SF.fit.FitSpec_Off = FitSpec;
Residual_FitOFF = real(OffSpecData(Freq))' - real(FitSpec)';
z1 = abs(JET_STR.SR.data.frequency_axis - Cho_lowerbound);
y1 = find(min(z1) == z1);
z2 = abs(JET_STR.SR.data.frequency_axis - Cho_upperbound);
y2 = find(min(z2) == z2);
JET_STR.SF.fit.Cho_GOF = ...
    sum(Residual_FitOFF(y1 - CrChoNAA_freqbounds(1) + 1 : y2 - CrChoNAA_freqbounds(1) + 1) .^ 2) / ...
    sum(real(OffSpecData(CrChoNAA_freqbounds(1) + y1 - y2 : CrChoNAA_freqbounds(1))) .^2 );

z1 = abs(JET_STR.SR.data.frequency_axis - Cr_lowerbound);
y1 = find(min(z1) == z1);
z2 = abs(JET_STR.SR.data.frequency_axis - Cr_upperbound);
y2 = find(min(z2) == z2);
JET_STR.SF.fit.Cr_GOF = ...
    sum(Residual_FitOFF(y1 - CrChoNAA_freqbounds(1) + 1 : y2 - CrChoNAA_freqbounds(1) + 1) .^2 ) / ...
    sum(real(OffSpecData(y1 + 1 : y2 + 1)) .^2 );
z1 = abs(JET_STR.SR.data.frequency_axis - NAA_lowerbound);
y1 = find(min(z1) == z1);
z2 = abs(JET_STR.SR.data.frequency_axis - NAA_upperbound);
y2 = find(min(z2) == z2);
JET_STR.SF.fit.NAA_GOF = ...
    sum(Residual_FitOFF(y1 - CrChoNAA_freqbounds(1) + 1 : y2 - CrChoNAA_freqbounds(1) + 1) .^2 ) / ...
    sum(real(OffSpecData(y2 : CrChoNAA_freqbounds(end))) .^2 );

h514 = figure(514);
subplot(2, 2, 3:4)
plot(JET_STR.SR.data.frequency_axis(DiffShow_freqbounds), real(OffSpecData(DiffShow_freqbounds)));
hold on;
plot(JET_STR.SR.data.frequency_axis(Freq), [real(FitSpec)', real(OffSpecData(Freq))' - real(FitSpec)']);
legend('Data','Fit','Residual');
axis([JET_STR.SF.params.Display_lowerbound, JET_STR.SF.params.Display_upperbound, ...
    min(0, 1 * min(real(OffSpecData(Freq)))), ...
    2 * max(real(OffSpecData(Freq)))]);
set(gca, 'XDir', 'reverse');
title('Cho Cr NAA Fitted');
set(gcf, 'color', 'w');

Cr_Amplitude = FitParams(1, 1); 
%Cr_R2a = FitParams(1, 2);
%Cr_R2b = FitParams(1, 3);
%Cr_ChemShift = FitParams(1, 4);
%if size(FitParams, 2) == 5
%    Cr_Phase = FitParams(2, 5);
%elseif size(FitParams, 2) == 4
%    Cr_Phase = 0;
%end

Cho_Amplitude = FitParams(2, 1);
%Cho_R2a = FitParams(2, 2);
%Cho_R2b = FitParams(2, 3);
%Cho_ChemShift = FitParams(2,4);
%if size(FitParams, 2) == 5
%    Cho_Phase = FitParams(2, 5);
%elseif size(FitParams,2) == 4
%    Cho_Phase = 0;
%end

NAA_Amplitude = FitParams(3, 1);
%NAA_R2a = FitParams(3, 2);
%NAA_R2b = FitParams(3, 3);
%NAA_ChemShift = FitParams(3, 4);
%if size(FitParams, 2) == 5
%    NAA_Phase = FitParams(3, 5);
%elseif size(FitParams, 2) == 4
%    NAA_Phase = 0;
%end

%baseline0 = FitParams(4, 1);
%baseline1 = FitParams(4, 2);


%% get the metabolite concentrations (the amplitude of FID)
GABA2Cr_Ratio = GABA_Amplitude / Cr_Amplitude;
Glx2Cr_Ratio = Glx_Amplitude / Cr_Amplitude;                                                                           

% save into the JET_STR.
JET_STR.SF.fit.GABA = GABA_Amplitude;
JET_STR.SF.fit.Glx = Glx_Amplitude;
JET_STR.SF.fit.Cr = Cr_Amplitude;
JET_STR.SF.fit.Cho = Cho_Amplitude;
JET_STR.SF.fit.NAA = NAA_Amplitude;

JET_STR.SF.fit.GABA2Cr_Ratio = GABA2Cr_Ratio;
JET_STR.SF.fit.Glx2Cr_Ratio = Glx2Cr_Ratio;

if isfield(JET_STR, 'HG_lowerbound')
    HG2Cr_Ratio = HG_Amplitude / Cr_Amplitude;
    JET_STR.SF.fit.HG = HG_Amplitude;
    JET_STR.SF.fit.HG2Cr_Ratio = HG2Cr_Ratio;
end

%%
Subject_foldername = JET_STR.Loader.subject_foldername_all(current_subject_folder_index).name;

%d_datetime = datetime('today');
%DateString = datestr(d_datetime);
DateString = datestr(now, 'yyyy-mm-dd');
%Report_dir = strcat('Vasile_Report/JET_Report_',DateString, '/', Subject_foldername);
Report_dir = strcat(JET_STR.Report.report_prefix, DateString, '/', Subject_foldername);
%mkdir(Report_dir)

if (exist(Report_dir, 'dir') ~= 7)
    eval(['mkdir ', Report_dir]);
end
saveas(h514, strcat(Report_dir, '/', Subject_foldername, '_Fitting_plot.png'))

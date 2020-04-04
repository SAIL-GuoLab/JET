function JET_STR = JET_Report(JET_STR, current_subject_folder_index)
%JET_REPORT Last step in the JET pipeline: generate the report.
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

close all

OffSpecData = fftshift(fft(ifft(fftshift(JET_STR.SR.data.off_average_SR)) .* ...
    (exp( - (JET_STR.SR.data.time_zeropad) * JET_STR.SF.params.LB * pi))));
DiffSpecData = fftshift(fft(ifft(fftshift(JET_STR.SR.data.diff_average_SR)) .* ...
    (exp( - (JET_STR.SR.data.time_zeropad) * JET_STR.SF.params.LB * pi))));

h200 = figure(200);

SpecDisplay_freqbounds = ...
JET_helper_function_find_freqbounds(JET_STR.SR.data.frequency_axis,...
4.2, 1.5);

set(gcf, 'units', 'normalized', 'outerposition', [0, 0, 0.5, 1]);
set(0, 'DefaulttextInterpreter', 'none')
gabafile = JET_STR.Loader.MEGAPRESS_filelist{1};
subjectID = gabafile(1 : end - 15);
subjectName_cell = strsplit(subjectID, '/');
Subject_foldername = subjectName_cell{end};

subplot(6, 6, [1, 2, 7, 8])
imagesc(real(JET_STR.SR.data.on_all_more_processing_CombChannel_SepRep_SR_final(:, SpecDisplay_freqbounds))');
% imagesc(abs(JET_STR.SR.data.on_all(SpecDisplay_freqbounds,:)));
xlabel('Repetition')
ylabel('ppm');
yticks(1 : round(length(SpecDisplay_freqbounds) / 5) : 5 * round(length(SpecDisplay_freqbounds) / 5))
yticklabels(num2cell(round(JET_STR.SR.data.frequency_axis(SpecDisplay_freqbounds(1 : round(length(SpecDisplay_freqbounds) / 5) : 5 * round(length(SpecDisplay_freqbounds) / 5))), 2)));
set(gca, 'fontsize', 8); ytickangle(90)
axis square
title('On');

subplot(6, 6, [3, 4, 9, 10])
imagesc(real(JET_STR.SR.data.on_all_more_processing_CombChannel_SepRep_SR_final(:, SpecDisplay_freqbounds))')
% imagesc(abs(JET_STR.SR.data.off_all(SpecDisplay_freqbounds,:)));
xlabel('Repetition')
ylabel('ppm');
yticks(1 : round(length(SpecDisplay_freqbounds) / 5) : 5 * round(length(SpecDisplay_freqbounds) / 5))
yticklabels(num2cell(round(JET_STR.SR.data.frequency_axis(SpecDisplay_freqbounds(1 : round(length(SpecDisplay_freqbounds) / 5) : 5 * round(length(SpecDisplay_freqbounds) / 5))), 2)));
set(gca, 'fontsize', 8); ytickangle(90)
axis square
title('Off')

subplot(6, 6, [5, 6, 11, 12])
imagesc(real(JET_STR.SR.data.diff_all_more_processing_CombChannel_SepRep_SR_final(:, SpecDisplay_freqbounds))')
% imagesc(abs(JET_STR.SR.data.diff_all(SpecDisplay_freqbounds,:)));
xlabel('Repetition')
ylabel('ppm');
yticks(1 : round(length(SpecDisplay_freqbounds) / 5) : 5 * round(length(SpecDisplay_freqbounds) / 5))
yticklabels(num2cell(round(JET_STR.SR.data.frequency_axis(SpecDisplay_freqbounds(1 : round(length(SpecDisplay_freqbounds) / 5) : 5 * round(length(SpecDisplay_freqbounds) / 5))), 2)));
set(gca,'fontsize',8);ytickangle(90);
axis square
title('Diff')

subplot(6, 6, 13 : 18)
plot(JET_STR.SR.data.icoshift_drift_Hz(:, 1), 'b-V')
hold on
plot(JET_STR.SR.data.icoshift_drift_Hz(:, 2), 'g-V')
legend('On','Off')
ylabel('Frequency Drift (Hz)')
% xlabel('Repetition')

subplot(6, 6, 19 : 24)
plot(JET_STR.SR.data.ACME_drift_rad(:, 1), 'b-V')
hold on
plot(JET_STR.SR.data.ACME_drift_rad(:, 2), 'g-V')
legend('On','Off')
ylabel('Phase Drift (rad)')
% xlabel('Repetition')

subplot(6,6,25:30)
plot(JET_STR.Report.SNR.OnSpectra_SNR_DB, 'b-*')
hold on
plot(JET_STR.Report.SNR.OffSpectra_SNR_DB, 'g-*')
legend('On', 'Off')
ylabel('SNR (dB)')
% xlabel('Repetition')

subplot(6, 6, 31 : 36)
plot(JET_STR.SR.data.on_all_more_processing_SR_RepSmooth_similarity, 'b-.')
hold on
plot(JET_STR.SR.data.off_all_more_processing_SR_RepSmooth_similarity, 'g-.')
legend('On', 'Off')
ylabel('R to Template')
xlabel('Repetition')

set(gcf, 'color', 'w')

%% save the figures as report
%d_datetime = datetime('today');
%DateString = datestr(d_datetime);

%saveas(h200,strcat(JET_STR.Report.report_dir,'/',Subject_foldername,'_Loading.png'))
saveas(h200, strcat(JET_STR.Report.report_dir, '/', Subject_foldername, '_Loading.png'))

%% fitting report
h201 = figure(201);
set(gcf, 'units', 'normalized', 'outerposition', [0, 0, 0.5, 1]);

% Compute the frequency bounds.
CrChoNAA_freqbounds = ...
JET_helper_function_find_freqbounds(JET_STR.SR.data.frequency_axis,...
JET_STR.SR.params.CrChoNAA_upperbound, JET_STR.SR.params.CrChoNAA_lowerbound);

GABA_freqbounds = ...
JET_helper_function_find_freqbounds(JET_STR.SR.data.frequency_axis,...
JET_STR.SR.params.GABA_upperbound, JET_STR.SR.params.GABA_lowerbound);

Glx_freqbounds = ...
JET_helper_function_find_freqbounds(JET_STR.SR.data.frequency_axis,...
JET_STR.SR.params.Glx_upperbound, JET_STR.SR.params.Glx_lowerbound);

Display_upperbound = 6; Display_lowerbound = 0;
Display_freqbounds = ...
JET_helper_function_find_freqbounds(JET_STR.SR.data.frequency_axis, Display_upperbound, Display_lowerbound);

subplot(6, 6, [1, 2, 3, 7, 8, 9]);
plot(JET_STR.SR.data.frequency_axis(:), real(OffSpecData(:)), 'k')
hold on
plot(JET_STR.SR.data.frequency_axis(CrChoNAA_freqbounds), JET_STR.SF.fit.FitSpec_Off, 'm');
% plot(JET_STR.SR.data.frequency_axis(CrChoNAA_freqbounds),real(OffSpecData(CrChoNAA_freqbounds))-JET_STR.SF.fit.FitSpec_Off,'k');
legend('Off', 'On', 'Off Fit')
axis([Display_lowerbound, Display_upperbound,...
    -0.5 * max(real(OffSpecData(CrChoNAA_freqbounds))), ...
    1.5 * max(real(OffSpecData(CrChoNAA_freqbounds)))])
set(gca, 'XDir', 'reverse')
ylabel('amplitude (a.u.)')
title('Off Spectra Fitting Result')
set(gcf, 'color', 'w')

subplot(6, 6, [4, 5, 6, 10, 11, 12]);
plot(JET_STR.SR.data.frequency_axis(:), real(DiffSpecData(:)), 'b')
hold on
plot(JET_STR.SR.data.frequency_axis(Glx_freqbounds), JET_STR.SF.fit.FitSpec_Diff_Glx, 'm');
hold on
plot(JET_STR.SR.data.frequency_axis(GABA_freqbounds), JET_STR.SF.fit.FitSpec_Diff_GABA, 'r');
% plot(JET_STR.SR.data.frequency_axis(GABAGlx_freqbounds),real(DiffSpecData(GABAGlx_freqbounds))-JET_STR.SF.fit.FitSpec_Diff,'k');
legend('Data', 'Glx Fit', 'Glx Fit', 'Location', 'southwest')
axis([Display_lowerbound, Display_upperbound, ...
    1.1 * min(real(DiffSpecData(Display_freqbounds))), ...
    1.5 * max(real(DiffSpecData(Display_freqbounds)))])
set(gca, 'XDir', 'reverse')
% ylabel('amplitude (a.u.)')
title('Diff Spectra Fitting Result')
set(gcf, 'color', 'w')

subplot(6, 6, [13, 14, 15]);
plot(JET_STR.SR.data.frequency_axis(CrChoNAA_freqbounds), ...
    real(OffSpecData(CrChoNAA_freqbounds)) - JET_STR.SF.fit.FitSpec_Off, 'k');
hold on;
plot(linspace(Display_lowerbound, Display_upperbound), ...
    0 .* linspace(Display_lowerbound,Display_upperbound), 'r-.')
axis([Display_lowerbound, Display_upperbound, ...
    1.5 .* min(real(OffSpecData(CrChoNAA_freqbounds)) - JET_STR.SF.fit.FitSpec_Off) ...
    1.5 .* max(real(OffSpecData(CrChoNAA_freqbounds)) - JET_STR.SF.fit.FitSpec_Off)])
set(gca, 'XDir', 'reverse')
legend('Residual', 'Location', 'southwest')
xlabel('ppm')
ylabel('amplitude (a.u.)')

subplot(6, 6, [16, 17, 18]);
plot(JET_STR.SR.data.frequency_axis(Glx_freqbounds), ...
    real(DiffSpecData(Glx_freqbounds)) - JET_STR.SF.fit.FitSpec_Diff_Glx, 'k');
hold on;
plot(JET_STR.SR.data.frequency_axis(GABA_freqbounds), ...
    real(DiffSpecData(GABA_freqbounds)) - JET_STR.SF.fit.FitSpec_Diff_GABA, 'k');
hold on;
plot(linspace(Display_lowerbound, Display_upperbound), ...
    0 .* linspace(Display_lowerbound, Display_upperbound), 'r-.')
axis([Display_lowerbound, Display_upperbound, ...
    1.5 .* min(real(DiffSpecData(Glx_freqbounds)) - JET_STR.SF.fit.FitSpec_Diff_Glx) ...
    1.5 .* max(real(DiffSpecData(Glx_freqbounds)) - JET_STR.SF.fit.FitSpec_Diff_Glx)])
set(gca, 'XDir', 'reverse')
legend('Residual', 'Location', 'southwest')
xlabel('ppm')
% ylabel('amplitude (a.u.)')

% Off Part

subplot(6, 6, 19 : 24);
text(0.05, 0, 'NAA');
text(0.2, 0, num2str(JET_STR.SF.fit.NAA));
text(0.35, 0, num2str(JET_STR.Report.SNR.NAA_SNR_DB)); axis off
% text(0.5,0,num2str(JET_STR.SF.fit.NAA_PeakHeight_SNR));axis off
% text(0.65,0,num2str(JET_STR.SF.fit.NAA_Area_SNR));axis off
text(0.8,0,strcat(num2str(100*JET_STR.SF.fit.NAA_GOF),'%')); axis off

text(0.05, 0.2, 'Cr');
text(0.2, 0.2, num2str(JET_STR.SF.fit.Cr));
text(0.35, 0.2, num2str(JET_STR.Report.SNR.ChoCr_SNR_DB)); axis off
% text(0.5,0.2,num2str(JET_STR.SF.fit.Cr_PeakHeight_SNR));axis off
% text(0.65,0.2,num2str(JET_STR.SF.fit.Cr_Area_SNR));axis off
text(0.8,0.2,strcat(num2str(100*JET_STR.SF.fit.Cr_GOF),'%')); axis off

text(0.05, 0.4, 'Cho');
text(0.2, 0.4, num2str(JET_STR.SF.fit.Cho));
text(0.35, 0.4, num2str(JET_STR.Report.SNR.ChoCr_SNR_DB)); axis off
% text(0.5,0.4,num2str(JET_STR.SF.fit.Cho_PeakHeight_SNR));axis off
% text(0.65,0.4,num2str(JET_STR.SF.fit.Cho_Area_SNR));axis off
text(0.8, 0.4, strcat(num2str(100 * JET_STR.SF.fit.Cho_GOF), '%')); axis off

h = text(0.05, 0.6, 'Metabolite');
set(h, 'Color', [1, 0, 0])
h = text(0.2, 0.6, 'Amplitude');
set(h, 'Color', [0, 0, 1])
h = text(0.35, 0.6, 'SNR(dB)'); axis off
set(h, 'Color', [0, 0, 1])
h = text(0.5, 0.6, 'SNR(H/std(N))'); axis off
set(h, 'Color', [0, 0, 1])
h = text(0.65, 0.6, 'SNR(A/std(N))'); axis off
set(h, 'Color', [0, 0, 1])
h = text(0.8, 0.6, 'GOF(NMSE)'); axis off
set(h, 'Color', [0, 0, 1])

% Diff Part
subplot(6, 6, 25 : 30);
text(0.05, 0.2, 'GABA');
text(0.2, 0.2, num2str(JET_STR.SF.fit.GABA));
text(0.35, 0.2, num2str(JET_STR.Report.SNR.GABA_SNR_DB)); axis off
% text(0.5,0.2,num2str(JET_STR.SF.fit.GABA_PeakHeight_SNR));axis off
% text(0.65,0.2,num2str(JET_STR.SF.fit.GABA_Area_SNR));axis off
text(0.8, 0.2, strcat(num2str(100 * JET_STR.SF.fit.GABA_GOF), '%')); axis off

text(0.05, 0.4, 'Glx');
text(0.2, 0.4, num2str(JET_STR.SF.fit.Glx));
text(0.35, 0.4, num2str(JET_STR.Report.SNR.Glx_SNR_DB)); axis off
% text(0.5,0.4,num2str(JET_STR.SF.fit.Glx_PeakHeight_SNR));axis off
% text(0.65,0.4,num2str(JET_STR.SF.fit.Glx_Area_SNR));axis off
text(0.8, 0.4, strcat(num2str(100 * JET_STR.SF.fit.Glx_GOF), '%')); axis off

h = text(0.05, 0.6, 'Metabolite');
set(h, 'Color', [1, 0, 0])
h = text(0.2, 0.6, 'Amplitude');
set(h, 'Color', [0, 0, 1])
h = text(0.35, 0.6, 'SNR(dB)'); axis off
set(h, 'Color', [0, 0, 1])
h = text(0.5, 0.6, 'SNR(H/std(N))'); axis off
set(h, 'Color', [0, 0, 1])
h = text(0.65, 0.6, 'SNR(A/std(N))'); axis off
set(h, 'Color', [0, 0, 1])
h = text(0.8, 0.6, 'GOF(NMSE)'); axis off
set(h, 'Color', [0, 0, 1])

% Diff Part
subplot(6, 6, 31 : 36);
text(0.05, 0.2, 'GABA/Cho');
text(0.2, 0.2, num2str(JET_STR.SF.fit.GABA ./ JET_STR.SF.fit.Cho));
text(0.35, 0.2, 'GABA/Cr'); axis off
text(0.5, 0.2, num2str(JET_STR.SF.fit.GABA ./ JET_STR.SF.fit.Cr)); axis off
text(0.65, 0.2, 'GABA/NAA'); axis off
text(0.8, 0.2, num2str(JET_STR.SF.fit.GABA ./ JET_STR.SF.fit.NAA)); axis off

text(0.05, 0.4, 'Glx/Cho');
text(0.2, 0.4, num2str(JET_STR.SF.fit.Glx ./ JET_STR.SF.fit.Cho));
text(0.35, 0.4, 'Glx/Cr'); axis off
text(0.5, 0.4, num2str(JET_STR.SF.fit.Glx ./ JET_STR.SF.fit.Cr)); axis off
text(0.65, 0.4, 'Glx/NAA'); axis off
text(0.8, 0.4, num2str(JET_STR.SF.fit.Glx ./ JET_STR.SF.fit.NAA)); axis off

h = text(0.05, 0.6, 'Metabolite');
set(h, 'Color', [1, 0, 0])
h = text(0.2, 0.6, 'Amplitude Ratio');
set(h, 'Color', [0, 0, 1])
h = text(0.35, 0.6, 'Metabolite'); axis off
set(h, 'Color', [1, 0, 0])
h = text(0.5, 0.6, 'Amplitude Ratio'); axis off
set(h, 'Color', [0, 0, 1])
h = text(0.65, 0.6, 'Metabolite'); axis off
set(h, 'Color', [1, 0, 0])
h = text(0.8, 0.6, 'Amplitude Ratio'); axis off
set(h, 'Color', [0, 0, 1])

%% save the figures as report

%saveas(h201,strcat(JET_STR.Report.report_dir,'/',Subject_foldername,'_Fitting.png'))
saveas(h201, strcat(JET_STR.Report.report_dir, '/', Subject_foldername, '_Fitting.png'))

%% temp code

%Spec_off_real = real(OffSpecData);
%Spec_diff_real=real(DiffSpecData);

%Spec_off_real_all=real(JET_STR.SR.data.off_all);
%Spec_diff_real_all=real(JET_STR.SR.data.diff_all);

%Spec_off(current_subject_folder_index,:)=OffSpecData;
%Spec_diff(current_subject_folder_index,:)=DiffSpecData                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   ;

% get the FWHM of Cr pic
[Cr_freqbound, Cr_upperbound_location, Cr_lowerbound_location] = ...
JET_helper_function_find_freqbounds(JET_STR.SR.data.frequency_axis,...
JET_STR.SR.params.Cr_upperbound, JET_STR.SR.params.Cr_lowerbound + 0.2);

Cr_Amp = max(real(OffSpecData(Cr_freqbound)));
z = abs(real(OffSpecData) - Cr_Amp);
Cr_idx = find(min(z) == z);

z = abs(real(OffSpecData(Cr_lowerbound_location : Cr_idx) - Cr_Amp /2));
Deltamin = find(min(z) == z) + Cr_lowerbound_location;
z = abs(real(OffSpecData(Cr_idx : Cr_upperbound_location) - Cr_Amp / 2));
Deltamax = find(min(z)==z) + Cr_idx;

Deltappm = abs(JET_STR.SR.data.frequency_axis(Deltamax) - JET_STR.SR.data.frequency_axis(Deltamin));
DeltaHz_Cr = Deltappm * 400 ;

% get the metabolite amplitudes
GABA_Amplitude(current_subject_folder_index) = JET_STR.SF.fit.GABA;
Glx_Amplitude(current_subject_folder_index) = JET_STR.SF.fit.Glx;
Cr_Amplitude(current_subject_folder_index) = JET_STR.SF.fit.Cr;
Cho_Amplitude(current_subject_folder_index) = JET_STR.SF.fit.Cho;
NAA_Amplitude(current_subject_folder_index) = JET_STR.SF.fit.NAA;
GABA_GOF(current_subject_folder_index) = JET_STR.SF.fit.GABA_GOF;
Glx_GOF(current_subject_folder_index) = JET_STR.SF.fit.Glx_GOF;
Cr_GOF(current_subject_folder_index) = JET_STR.SF.fit.Cr_GOF;
Cho_GOF(current_subject_folder_index) = JET_STR.SF.fit.Cho_GOF;
NAA_GOF(current_subject_folder_index) = JET_STR.SF.fit.NAA_GOF;
Cr_FWHM(current_subject_folder_index) = DeltaHz_Cr;
NSDR_GABA(current_subject_folder_index) = JET_STR.SF.fit.NSDR_GABA / GABA_Amplitude(current_subject_folder_index);
NSDR_Glx(current_subject_folder_index) = JET_STR.SF.fit.NSDR_Glx / Glx_Amplitude(current_subject_folder_index);

%%
subject_name = JET_STR.Loader.subject_foldername_all(current_subject_folder_index).name;
table =[string(subject_name), string(GABA_Amplitude(current_subject_folder_index)), ...
    string(Glx_Amplitude(current_subject_folder_index)), string(Cr_Amplitude(current_subject_folder_index)), ...
    string(Cho_Amplitude(current_subject_folder_index)), string(NAA_Amplitude(current_subject_folder_index)), ...
    string(GABA_GOF(current_subject_folder_index)), string(Glx_GOF(current_subject_folder_index)), ...
    string(Cr_GOF(current_subject_folder_index)), string(Cho_GOF(current_subject_folder_index)), ...
    string(NAA_GOF(current_subject_folder_index)), string(Cr_FWHM(current_subject_folder_index)), ...
    string(NSDR_GABA(current_subject_folder_index)), string(NSDR_Glx(current_subject_folder_index))];
file = fopen(JET_STR.Report.table.filename, 'a');
fprintf(file, '\n');
fprintf(file, '%s,', table);
fclose(file);
end
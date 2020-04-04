function JET_STR = JET_Loader_preinitialize(JET_STR)
%JET_Loader_preinitialize A helper function to facilitate JET_Loader, which
%is the first step of the JET pipeline. This will set some initial
%parameters, but if any of these values is detected in the header when we
%load the actual file, it will be overwritten.
%
% Input arguments
% - JET_STR : The JET struct for the study.
%
% Output arguments
% - JET_STR : The modified JET struct for the study.
%
% Chen "Raphael" Liu (cl3760@columbia.edu) & Jia Guo (jg3400@columbia.edu),
% 02/24/2020.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Spectrum Registration Parameters %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    JET_STR.SR.params.SR_iterations = 2; % Iterations of spectrum registration.
    JET_STR.SR.params.field_strength = 3; % 3 Tesla.
    JET_STR.SR.params.water_ppm = 4.7; % Water resonates at 4.7 ppm.
    JET_STR.SR.params.LarmorFreq = 400; % Equals to 42.5781 * field_strength, in our case equals to 42.5781*9.4 = 400. This should be parsed from headers where possible.
    JET_STR.SR.params.target = 'GABA'; % Other options are 'GSH' and 'Glx' (now also implemented).
    JET_STR.SR.params.ONOFForder = 'onfirst'; %Options are 'onfirst' or 'offfirst'.
    %JET_STR.SR.params.SegMode = 1; % 1 = average the segment then process the data, 2 = directly process the data.
    JET_STR.SR.params.LS_mode = 1;
    
    JET_STR.SR.params.Freq_Shift_Th = 13;
    
    %JET_STR.SR.params.LB4FC = 12;
    %JET_STR.SR.params.LB4SR_additional = 12;
    JET_STR.SR.params.LB = 12;
    JET_STR.SR.params.LB_water = 0;
    JET_STR.SR.params.ZeroFillTo = 2^15;
    
    %AlignTo planned options: Cr; Cho; NAA; H20; CrOFF
    %JET_STR.Loader.params.AlignTo = 'Bruker'; % SpecReg default and recommended.

    % Options for initial Spectrum Registration using ACME and icoshift.
    JET_STR.SR.params.icoshift_mode = 'average'; % 'average' 'median' 'max'
    JET_STR.SR.params.icoshift_spectra = 'Off'; % 'Sum' 'Off'

    JET_STR.SR.params.icoshift_lowerbound1 = -10;
    JET_STR.SR.params.icoshift_upperbound1 = 10;

    JET_STR.SR.params.icoshift_lowerbound2 = 7;
    JET_STR.SR.params.icoshift_upperbound2 = 10;

    JET_STR.SR.params.SR_lowerbound = 0;
    JET_STR.SR.params.SR_upperbound = 7.5;

    JET_STR.SR.params.SR_RepSmooth_Alpha = 10;

    JET_STR.SR.params.SR_with_real_or_complex = 'complex';

    JET_STR.SR.params.SRfreq_lowerbound_On_to_Off_1 = -1;
    JET_STR.SR.params.SRfreq_upperbound_On_to_Off_1 = 0.2;

    JET_STR.SR.params.SRfreq_lowerbound_On_to_Off_2 = 2.7;
    JET_STR.SR.params.SRfreq_upperbound_On_to_Off_2 = 2.9;    
    
    JET_STR.SR.params.SRfreq_lowerbound_On_to_Off_3 = 3.2;
    JET_STR.SR.params.SRfreq_upperbound_On_to_Off_3 = 3.5;
    
    JET_STR.SR.params.SRfreq_lowerbound_On_to_Off_4 = 4.5;
    JET_STR.SR.params.SRfreq_upperbound_On_to_Off_4 = 10;
        
    JET_STR.SR.params.SRSpec_water_lowerbound = 4.75-3;
    JET_STR.SR.params.SRSpec_water_upperbound = 4.75+3;

    JET_STR.SR.params.Spectra_lowerbound = 1;
    JET_STR.SR.params.Spectra_upperbound = 4;

    JET_STR.SR.params.Display_lowerbound = -1;
    JET_STR.SR.params.Display_upperbound = 7;
    
    JET_STR.SR.params.Noise_lowerbound = 6.5;
    JET_STR.SR.params.Noise_upperbound = 7.8;
    
    JET_STR.SR.params.PhC_ACME_lowerbound = 1.5;
    JET_STR.SR.params.PhC_ACME_upperbound = 4;
    
    JET_STR.SR.params.PhC_ACMH_lowerbound = 3.5;
    JET_STR.SR.params.PhC_ACMH_upperbound = 4;
    
    JET_STR.SR.params.icoshift_water_lowerbound = 4.7-2;
    JET_STR.SR.params.icoshift_water_upperbound = 4.7+2;
    
    JET_STR.SR.params.PhC_water_lowerbound = 4.7-2;
    JET_STR.SR.params.PhC_water_upperbound = 4.7+2;
    
    % Options for fine-tuning Spectrum Registration.
    JET_STR.SR.params.water_lowerbound = 4.7-2;
    JET_STR.SR.params.water_upperbound = 4.7+2;
    
    JET_STR.SR.params.GABA_lowerbound = 2.7;
    JET_STR.SR.params.GABA_upperbound = 3.3;
    
    JET_STR.SR.params.Glx_lowerbound = 3.65;
    JET_STR.SR.params.Glx_upperbound = 3.85;
    
    JET_STR.SR.params.GlxGABA_lowerbound = 2.5;
    JET_STR.SR.params.GlxGABA_upperbound = 4.5;
    
    JET_STR.SR.params.NAA_lowerbound = 1.6;
    JET_STR.SR.params.NAA_upperbound = 2.2;
    JET_STR.SR.params.NAA_target_lowerbound = 1.5;
    JET_STR.SR.params.NAA_target_upperbound = 2.5;

    JET_STR.SR.params.Cr_lowerbound = 2.6;
    JET_STR.SR.params.Cr_upperbound = 3.11;

    JET_STR.SR.params.ChoCr_lowerbound = 2.8;
    JET_STR.SR.params.ChoCr_upperbound = 3.5;

    JET_STR.SR.params.Cho_lowerbound = 3.1;
    JET_STR.SR.params.Cho_upperbound = 3.6;
    
    JET_STR.SR.params.CrChoNAA_lowerbound = 1.5;
    JET_STR.SR.params.CrChoNAA_upperbound = 3.6;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Spectral Fitting Parameters %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    JET_STR.SF.params.LB = 8;
    JET_STR.SF.params.LB_metabolites = 2;
    JET_STR.SF.params.ZeroFillTo = 2^15;

    JET_STR.SF.params.Display_lowerbound = 1.5;
    JET_STR.SF.params.Display_upperbound = 4.5;
    
    JET_STR.SF.params.GABA_lowerbound = 2.7;
    JET_STR.SF.params.GABA_upperbound = 3.3;
    
    JET_STR.SF.params.Glx_lowerbound = 3.65;
    JET_STR.SF.params.Glx_upperbound = 3.85;
    
    JET_STR.SF.params.GlxGABA_lowerbound = 2.5;
    JET_STR.SF.params.GlxGABA_upperbound = 4.5;

    JET_STR.SF.params.Cr_lowerbound = 1.5;
    JET_STR.SF.params.Cr_upperbound = 3.1;
    
    JET_STR.SF.params.Cho_lowerbound = 3.1;
    JET_STR.SF.params.Cho_upperbound = 3.6;

    JET_STR.SF.params.NAA_lowerbound = 1.5;
    JET_STR.SF.params.NAA_upperbound = 2.4;
    
    JET_STR.SF.params.CrChoNAA_lowerbound = 1.5;
    JET_STR.SF.params.CrChoNAA_upperbound = 3.6;

    JET_STR.SF.params.extra_lowerbound1 = -1;
    JET_STR.SF.params.extra_upperbound1 = -0.5;
    
    JET_STR.SF.params.extra_lowerbound2 = 9;
    JET_STR.SF.params.extra_upperbound2 = 9.5;
    
    % Output Parameters
    JET_STR.Report.save_JET_STR.mat = 1; % 1 = YES, save JET_STR as .mat file.
    JET_STR.Report.save_JET_STR.sdat = 0; % 1 = YES, save JET_STR as .sdat file.
end

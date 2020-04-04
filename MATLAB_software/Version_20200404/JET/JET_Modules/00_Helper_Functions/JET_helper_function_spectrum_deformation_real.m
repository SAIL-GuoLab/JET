% These are some helper functions that are used in the JET pipeline.
% Some processing are used over and over again in multiple scripts or
% stages, such that we found it reasonable to allocate designated helper
% functions for them.
%
% Chen "Raphael" Liu, Vasile Stupar and Jia Guo, 02/24/2020.

function Spec_Deformed = JET_helper_function_spectrum_deformation_real(x, time, Freq, FID_Orig, LS_mode)
% Spectrum quantification with simulated basis set.
% Jia, 10/10/16
% Please contact jg3400@columbia.edu if you have questions.
% x(1,1) --- Amplitude
% x(1,2) --- R2
% x(1,3) --- ChemShift
% x(1,4) --- zero order Phase
% x(2,1) --- baseline

switch LS_mode
    case 1 % Lorenzian Line Broadening
        if size(x,2)==2 % 
            ChemShift = x(1,1);
            Phase = x(1,2);
            FID_Deformed =  1.*FID_Orig.*exp(-0.*time).*exp(1i.*time.*ChemShift.*300.*6.28+1i.*Phase);
        elseif size(x,2)==3 % 
            Amplitude = x(1,1);
            ChemShift = x(1,2);
            Phase = x(1,3);
            FID_Deformed = Amplitude.*FID_Orig.*exp(-0.*time).*exp(1i.*time.*ChemShift.*300.*6.28+1i.*Phase);
        elseif size(x,2)==4 % zero order phase
            Amplitude = x(1,1);
            R2 = x(1,2);
            ChemShift = x(1,3);
            Phase = x(1,4);
            FID_Deformed =  Amplitude.*FID_Orig.*exp(-R2.*time).*exp(1i.*time.*ChemShift.*300.*6.28+1i.*Phase);
        end
    case 2 % Gaussian Line Broadening
        if size(x,2)==2 % zero order phase
            ChemShift = x(1,1);
            Phase = x(1,2);
            FID_Deformed =  1.*FID_Orig.*exp(-0.*time.^2).*exp(1i.*time.*ChemShift.*300.*6.28+1i.*Phase);
        elseif size(x,2)==3 % 
            Amplitude = x(1,1);
            ChemShift = x(1,2);
            Phase = x(1,3);
            FID_Deformed = Amplitude.*FID_Orig.*exp(-0.*time).*exp(1i.*time.*ChemShift.*300.*6.28+1i.*Phase);
        elseif size(x,2)==4 % zero order phase
            Amplitude = x(1,1);
            R2 = x(1,2);
            ChemShift = x(1,3);
            Phase = x(1,4);
            FID_Deformed = Amplitude.*FID_Orig.*exp(-(R2.*time).^2).*exp(1i.*time.*ChemShift.*300.*6.28+1i.*Phase);
        end
end

Spec_Deformed_complex = fliplr(fft(FID_Deformed));
Spec_Deformed = real(Spec_Deformed_complex(Freq));

if size(x,1)==2
    Spec_Deformed = Spec_Deformed+x(end,1).*Freq+x(end,2);
end
end
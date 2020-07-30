function FitSpec = SQSBS(x, time, Freq, FID_basis, LS_mode, FrqRef) % spectrum quantification with simulated basis set

% Jia, 10/10/16
% Please contact jg3400@columbia.edu if you have questions.

switch LS_mode
    case 1 % Lorenzian Line Broadening
        if size(x,2)==3
            FitFID=0;
            for i=1:size(FID_basis,1)
                Amplitude = x(i,1);
                R2 = x(i,2);
                ChemShift = x(i,3);
                FitFID = FitFID + Amplitude.*FID_basis(i,:).*exp(-R2.*time).*exp(1i.*time.*ChemShift.*FrqRef.*2*pi);
            end
        elseif size(x,2)==4
            FitFID=0;
            for i=1:size(FID_basis,1)
                Amplitude = x(i,1);
                R2 = x(i,2);
                ChemShift = x(i,3);
                Phase = x(i,4);
                FitFID = FitFID + Amplitude.*FID_basis(i,:).*exp(-R2.*time).*exp(1i.*time.*ChemShift.*FrqRef.*2*pi+1i.*Phase);
            end
        end
    case 2 % Gaussian Line Broadening
        if size(x,2)==3
            FitFID=0;
            for i=1:size(FID_basis,1)
                Amplitude = x(i,1);
                R2 = x(i,2);
                ChemShift = x(i,3);
                FitFID = FitFID + Amplitude.*FID_basis(i,:).*exp(-(R2.*time).^2).*exp(1i.*time.*ChemShift.*FrqRef.*2*pi);
            end
        elseif size(x,2)==4
            FitFID=0;
            for i=1:size(FID_basis,1)
                Amplitude = x(i,1);
                R2 = x(i,2);
                ChemShift = x(i,3);
                Phase = x(i,4);
                FitFID = FitFID + Amplitude.*FID_basis(i,:).*exp(-(R2.*time).^2).*exp(1i.*time.*ChemShift.*FrqRef.*2*pi+1i.*Phase);
            end
        end
    case 3 % Voigt Line Broadening
        if size(x,2)==4
            FitFID=0;
            for i=1:size(FID_basis,1)
                Amplitude = x(i,1);
                R2a = x(i,2);
                R2b = x(i,3);
                ChemShift = x(i,4);
                FitFID = FitFID + Amplitude.*FID_basis(i,:).*exp(-R2a.*time).*exp(-(R2b.*time).^2).*exp(1i.*time.*ChemShift.*FrqRef.*2*pi);
            end
        elseif size(x,2)==5
            FitFID=0;
            for i=1:size(FID_basis,1)
                Amplitude = x(i,1);
                R2a = x(i,2);
                R2b = x(i,3);
                ChemShift = x(i,4);
                Phase = x(i,5);
                FitFID = FitFID + Amplitude.*FID_basis(i,:).*exp(-R2a.*time).*exp(-(R2b.*time).^2).*exp(1i.*time.*ChemShift.*FrqRef.*2*pi+1i.*Phase);
            end
        end
end

% FitSpec = fliplr(real(fftshift(fft(FitFID))));
FitSpec = real(fftshift(fft(FitFID))); % vs

% FitSpec = FitSpec(Freq);
FitSpec = FitSpec(Freq)+x(end,1).*Freq+x(end,2);
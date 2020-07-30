Spectral_ints=[1 100 100 140 140 200 200 250 250 300 300 550 550 670 670 800 800 1000 1000 1250 1250 1500 1500 1550 1550 1670 1670 1750 1750 1900 1900 2400 2400 2200 3000 3700];

load Mutant_waterL 
load Mutant_waterR 
load WT_waterL 
load WT_waterR

waterL=[Mutant_waterL;WT_waterL];

waterR=[Mutant_waterR;WT_waterR];

SpectralData=waterR;

[AlignedSpectralData_waterR, intervals, indexes] = icoshift ('average', SpectralData, Spectral_ints, 'f', [2 1 0], [1:1:size(SpectralData,2)]);

save AlignedSpectralData_waterR AlignedSpectralData_waterR
save Spectral_ints Spectral_ints



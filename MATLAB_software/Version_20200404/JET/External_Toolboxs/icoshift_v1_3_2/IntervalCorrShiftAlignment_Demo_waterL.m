%Spectral_ints=[1 32 33 169 170 359 360 499 500 649 650 749 750 909 910 980 980 1049 1050 1249 1250 1399 1400 1700 1700 2001];

Spectral_ints = [1 1500 1500 2100 2100 2300 2300 2430 2430 2500 2500 2550 2550 2600 2600 2900 2900 3500 3500 4200 4200 4600 4600 5301];

load Mutant_waterL 
load Mutant_waterR 
load WT_waterL 
load WT_waterR

waterL=[Mutant_waterL;WT_waterL];

waterR=[Mutant_waterR;WT_waterR];

SpectralData=waterL;

[AlignedSpectralData_waterL, intervals, indexes] = icoshift ('average', SpectralData, Spectral_ints, 'f', [2 1 0], [1:1:size(SpectralData,2)]);

save AlignedSpectralData_waterL AlignedSpectralData_waterL
save Spectral_ints Spectral_ints



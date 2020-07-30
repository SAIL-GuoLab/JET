%% Load demo Wine Data

load WineData
%% Plot raw data
% the region between 1.25 and 1.45 ppm contains the succinic acid signal

plot (ppm, X');
set(gca,'XDir','reverse');

%% Plot the lactic acid zoomed region (1.25 - 1.45 ppm)

plot (ppm(7151:7550), X(:,7151:7550));
axis tight
set(gca,'XDir','reverse');

%% iCOshift 1: aligns the whole spectra according to a reference signal selected (Ethanol CH3 resonance)
% Try to zoom on one subplot to view more features. Turn off the zoom
% button and try to click on one spectrum or on one interval to see more details.

[AlignedWineData, intervals, indexes] = icoshift ('average', X, (7551:7750), 'f',[2 1 0]);

%% iCOshift 2: splits the dataset in 50 regular intervals and aligns each of them separately
% COMMAND: [xCS,ints,ind] = icoshift(xT,xP,inter,{n},{options})
% [AlignedWineData, intervals, indexes] = icoshift (target sample, Dataset, intervals mode definition, allowed shift mode , plot options);

[AlignedWineData, intervals, indexes] = icoshift ('average', X, 50, 'f', [2 1 0]);

%% iCOshift 3: splits the dataset in regular intervals 800 points wide and search for the "best" allowed shift for each of them separately

[AlignedWineData, intervals, indexes] = icoshift ('average', X, '800', 'b', [2 1 0]);

%% iCOshift 4: splits the dataset in pre-defined intervals (on the basis of the user's knowledge) and aligns each of them
% The ppm scale is used for the final plot and the x axes is automatically reversed
% The vector used for defining the intervals is compatible with the
% algorithms in the iToolbox (e.g. iPLS and iPCA)
 
[AlignedWineData, intervals, indexes] = icoshift ('average', X, wine_ints, 'f', [2 1 0], ppm);

%% iCOshift 4a: Like the previous example but using an intermediate "coshift" step

[AlignedWineData, intervals, indexes] = icoshift ('max', X, wine_ints, 'f', [2 1 1], ppm);
set(gca,'Xlim',[1.28 1.44], 'Ylim', [-100 9e06]);


%% iCOshift 4b: Like the previous example but using also missing values (NaN) in the last setp

[AlignedWineData, intervals, indexes] = icoshift ('max', X, wine_ints, 'f', [2 0 1], ppm);
set(gca,'Xlim',[1.28 1.44], 'Ylim', [-100 9e06]);

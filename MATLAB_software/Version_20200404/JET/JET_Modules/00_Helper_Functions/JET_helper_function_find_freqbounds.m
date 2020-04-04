% These are some helper functions that are used in the JET pipeline.
% Some processing are used over and over again in multiple scripts or
% stages, such that we found it reasonable to allocate designated helper
% functions for them.
%
% Chen "Raphael" Liu, Vasile Stupar and Jia Guo, 02/24/2020.

function [freqbounds, upperbound_location, lowerbound_location] = JET_helper_function_find_freqbounds(frequency_spectrum, defined_upperbound, defined_lowerbound)
z = abs(frequency_spectrum - defined_upperbound);
upperbound_location = find(min(z) == z);
z = abs(frequency_spectrum - defined_lowerbound);
lowerbound_location = find(min(z) == z);
freqbounds = lowerbound_location : upperbound_location;
end
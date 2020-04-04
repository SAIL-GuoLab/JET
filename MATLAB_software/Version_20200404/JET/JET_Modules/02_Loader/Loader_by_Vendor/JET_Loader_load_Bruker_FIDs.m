function organized_FIDs = JET_Loader_load_Bruker_FIDs(gabafile, JET_STR)
% JET_LOADER_LOAD_BRUKER_FIDS A helper function that reads the FIDs.
%
% The objectives of this functions include:
% 1) read the FIDs
% 2) correctly organize the FIDs into On-Off interleave FIDs
% 3) store the resulting organized FIDs in JET_STR
%
% Chen "Raphael" Liu (cl3760@columbia.edu) and Jia Guo (jg3400@columbia.edu), 02/24/2020.

switch JET_STR.Loader.coil_type
    case 'multi_array_surface_coil'
        organized_FIDs = JET_Loader_load_Bruker_FIDs_MultiArraySurfaceCoil(gabafile, JET_STR);
    case 'volume_coil'
        organized_FIDs = JET_Loader_load_Bruker_FIDs_VolumeCoil(gabafile, JET_STR);
end

end

function sp_dp = dephase_fun(sp_ft,phc)
% function sp_dp = dephase_fun(sp_ft,phc0,phc1)
% written by Chen Li 
% input:   sp_ft - complex data after fourier transform
%         phc0  - the zero order phase correction
%         phc1  - the first order phase correction
% output:  sp_dp - spectral data after phase correction
%
% Modified by Mike Tyszka, Caltech BIC
size_phc=length(phc);
if(size_phc == 2)
    phc0=phc(1);
    phc1=phc(2);
    
    phc0 = phc0 * pi/180;              % convert degree to radian
    phc1 = phc1 * pi/180;              % convert degree to radian
    
    % m complex spectra (rows) with n samples (cols)
    [m,n]=size(sp_ft);
    
    % Normalized index vector
    a_num = linspace(0,1,n);
    
    % Calculate a(i,j)
    a = phc0 + a_num * phc1;
    a = repmat(a,[m 1]);
    
    % Apply phase function
    sp_dp = sp_ft .* exp(1i * a);

elseif (size_phc == 1)
    phc0=phc(1);
    phc0 = phc0 * pi/180;              % convert degree to radian
%     phc1 = phc1 * pi/180;              % convert degree to radian
    
    % m complex spectra (rows) with n samples (cols)
    [m,n]=size(sp_ft);
%     
%     % Normalized index vector
%     a_num = linspace(0,1,n);
    
    % Calculate a(i,j)
    a = phc0;
    a = repmat(a,[m 1]);
    
    % Apply phase function
    sp_dp = sp_ft .* exp(1i * a);

end
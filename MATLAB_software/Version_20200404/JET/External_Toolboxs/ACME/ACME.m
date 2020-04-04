    function [NMRsp_dp,phc]=ACME(NMRsp,phcf)
% Chen Li's automatic spectrum phasing function
%
% [NMRsp_dp,phc0,phc1]=ACME(NMRsp,phcf)
%
% Original function by Chen Li 
%
% input   NMRsp     - NMR spectral data (complex row vector)
%         phcf      - initial estimates of phc0 and phc1, usually phcf is set to [0 0] 
%
% output  NMRsp_dp  - dephased NMR spectral data (complex row vector)
%         phc0      - the zero order phase correction (degree)
%         phc1      - the first order phase correction(degree)
%
% The first derivative is used for calculation of entropy.
% subprograms: dephase_fun.m and entropy_fun.m
%
% REQUIRED REFERENCE: An efficient algorithm for automatic phase correction
% of NMR spectra based on entropy minimization", Li Chen, Zhiqiang Weng,
% LaiYoong Goh, and Marc Garland Journal of Magnetic Resonance 158(2002)
% 164-168

if nargin < 2 phcf = [0 0]; end

% Force row vector
NMRsp = NMRsp(:)';

fac = fminsearch('My_entropy_fun',phcf,[],NMRsp);

if length(fac)==1
    NMRsp_dp = dephase_fun(NMRsp,fac); 
    NMRsp_dp = conj(NMRsp_dp);
    phc = fac;
elseif length(fac)==2
    NMRsp_dp = dephase_fun(NMRsp,fac(1),fac(2));
    NMRsp_dp = conj(NMRsp_dp);
    phc(1) = fac(1);
    phc(2) = fac(2);
end
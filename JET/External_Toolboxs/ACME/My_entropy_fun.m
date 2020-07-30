function f = My_entropy_fun(x,s)
% function f = entropy_fun(x,s)
% written by Chen Li
% input:   x - PHC0 and PHC1
%          s - NMR Spectral data
% output   f - entropy value (Using the first derivative)
%
% Modifications: Mike Tyszka, Ph.D, Caltech BIC

% Derivative order for entropy calculation
dn = 2;

% Dephase
if length(x)==2
    phc0 = x(1);
    phc1 = x(2);
    s0 = dephase_fun(s,phc0,phc1);
elseif length(x)==1
    phc0 = x(1);
    s0 = dephase_fun(s,phc0);
end

s = real(s0);
si = imag(s0);
s_mag=(s.^2+si.^2).^0.5;

% Absolute Nth order derivative
ds_n = abs(diff(s,dn));
dsi_n = abs(diff(si,dn));

% Normalized derivative
p_n = ds_n ./ sum(ds_n);
p_ni = dsi_n ./ sum(dsi_n);

% Eliminate zeros (ln(1) = 0)
p_n(p_n == 0) = 1;
p_ni(p_ni == 0) = 1;

% Entropy measure
H_n = sum(-p_n .* log(p_n));
H_ni = sum(-p_ni .* log(p_ni));

% baseline approximation
baseline=linspace(s(1),s(end),length(s));
s=s-baseline;

% Calculation of penalty
Pfun = 0.0;
as = s - abs(s);
sumas = sum(sum(as));
if (sumas < 0)
    Pfun = Pfun + sum(sum((as./2).^2));
end

% area under the curve
area=sum(s);
areai=sum(si);
areaabs=sum(abs(s));
correlation=corr2(s,s_mag)*corr2(s,s_mag);
correlationi=corr2(abs(si),s_mag)*corr2(abs(si),s_mag);

% Penalty function to ensure positive bands
Penalty = 1000 * Pfun;
% Penalty = -correlation;

% % based on entrophy
% f = H_n;

f = H_n;

% % based on correlation to the abs spectrum and entrophy
% f = H_n + Penalty;

% % based on correlation to the abs spectrum
% f = -correlation;

% % based on area
f = H_n+abs(sum(si-min(abs(si))))/abs(sum(s-min(abs(s))));



function qvol = shelldens(bval, bprecision)
% qvol = shelldens(bval)
% 
%   Sampling density nonuniformity correction factor for multi-shell
%   samples, i.e. volume associated with each sample.
%
%   input:
%       bval: b-values, might have [-200, 200] deviation from prescribed values
%       bprecision: precesion to round up bval
%
%   output:
%       qvol: q-space sampling density correction factor
%
%   Reference:
%       [1] Tian Q, Yang G, Leuze C, Rokem A, Edlow BL, McNab JA, Generalized
%       Diffusion Spectrum Magnetic Resonance Imaging (GDSI) for Model-free
%       Reconstruction of the Ensemble Average Propagator. NeuroImage, 2019;
%       189: 497-515.
%
% (c) Qiyuan Tian, Stanford RSL, 2019 Jan

bval = bval(:); % column vector
bval_rnd = round(bval / bprecision) * bprecision; % round up b-value

bval_uniq = unique(bval_rnd); % unique b-value, double check 
disp('double check unique b-values are:');
disp(bval_uniq);
disp('otherwise change bprecision to correctly round up b-values.');

count = zeros(size(bval_uniq)); % number of each b-value
for ii = 1 : length(bval_uniq)
    count(ii) = sum(bval_rnd == bval_uniq(ii));
end

qval_uniq = sqrt(bval_uniq / max(bval_uniq)); % nomalized unique q-value, q~sqrt(b)

qval_contour = (qval_uniq(1 : end-1) + qval_uniq(2 : end)) / 2; % middle contour
qval_contour = [qval_contour; 2 * qval_uniq(end) - qval_contour(end)]; % extrapolate outer contour 

qvol_shell = diff(qval_contour .^ 3); % volume associated with each shell
qvol_shell = [qval_contour(1) .^ 3; qvol_shell]; % add in central sphere volume

qvol_shell = qvol_shell / qvol_shell(1); % normalize central sphere volume to 1
qvol_samp = qvol_shell ./ count; % volume associated with a sample on a shell

qvol = zeros(size(bval_rnd)); % qspace volume for each sample 
for ii = 1 : length(bval_uniq)
    b = bval_uniq(ii);
    qvol(bval_rnd == b) = qvol_samp(ii);
end



function [I_noisy_phase_stab, I_noisy_phase_stab_imag] = dwi_prepareFiles(file_magn,file_phase,opt)



fprintf('Preparing files ... \n')

% Load the files
I_noisy_magn = niftiread(file_magn);
if isempty(file_phase)
    fprintf('Magnitude-only data ...  \n')
    fprintf('Will use Rician VST framework ...  \n')
    I_noisy_phase = [];
else
    fprintf('Complex-valued data ... performing phase stabilization \n')
    I_noisy_phase = niftiread(file_phase);
end

% Some nii are saved as int format, so we convert to single or double
% depending on bits already there.
nii_info = niftiinfo(file_magn);
if nii_info.BitsPerPixel > 32
    I_noisy_magn = double(I_noisy_magn);
    I_noisy_phase = double(I_noisy_phase);
else
    I_noisy_magn = single(I_noisy_magn);
    I_noisy_phase = single(I_noisy_phase);
end

% Do the actual phase stabalization
if isempty(file_phase)
    I_noisy_phase_stab = I_noisy_magn;
    I_noisy_phase_stab_imag = [];
else
    I_noisy_phase_stab = perform_phase_stabilization(I_noisy_magn,I_noisy_phase,opt);
    clear I_noisy_magn I_noisy_phase
    I_noisy_phase_stab_imag = imag(I_noisy_phase_stab);
    I_noisy_phase_stab = real(I_noisy_phase_stab);
end

fprintf('Finished preparing files. \n\n');



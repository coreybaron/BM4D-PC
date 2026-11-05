function I_denoised = denoise_BM4D_PC(file_magn,file_phase,opt)
% ************************************************************************
%  This function implements the denoising method BM4D-PC[1]
%
%  % [1]. Vinicius P. Campos, Diego Szczupak, Tales Santini, Afonso C. Silva, Alessandro Foi, Marcelo A. C. Vieira, and Corey Baron.
%  % "BM4D-PC: nonlocal volumetric denoising of principal components of diffusion-weighted MR images"

% INPUTS: 
%
%            file_magn  : complete path to magnitude image nifti file (.nii or .nii.gz)
%
%                           OBS:  - The bval file must have .bval extension. 
%                                 - It must be in the same directory of the image to be denoised. 
%                                 - Its name must match magnitude image file name
%                                 - If no bval file is found, all dwi data will be considered for
%                                 noise estimation.
%         
%            file_phase : complete path to phase nifti image file (.nii or .nii.gz) - (optional)
%
%                        - if specified, phase stabilization will  be
%                        performed and denoise will occur in the complex domain
%                        
%                        - if not specified, only magnitude image will be
%                        used and denoising will be performed using the
%                        Rician VST framework.
%
%             opt        : struct containing optional parameters. Please check scripts - (optional)
%
%
% OUTPUT: 
%
%             I_denoised : the denoised image
%
%                            OBS: the algorithm already saves a denoised nifti file on the same folder of input file. 
%
%
% ****************** Usage examples ****************
%
% file_magn = 'path/to/magnitude_image.nii';
% file_phase = 'path/to/phase_image.nii';
%
% I_denoised = denoise_BM4D_PC(file_magn,file_phase);
%
% **************************************************
%
% file_magn = 'path/to/magnitude_image.nii';
%
% I_denoised = denoise_BM4D_PC(file_magn);
%
% **************************************************
%
% file_magn = 'path/to/magnitude_image.nii';
% opt.profile_bm4d = 'np';
%
% I_denoised = denoise_BM4D_PC(file_magn,[],opt);
%
% ************************************************************************
% AUTHOR:
%
%     Vinicius P. Campos
%     email: vpc24@pitt.edu
%
% ************************************************************************


if isempty(file_magn)
    error('Please provide magnitude file path')
end

if ~exist('opt','var') || ~exist('opt.profile_bm4d','var')
    opt.profile_bm4d = 'np';
end

if ~exist('file_phase','var')
    opt.flag_magnitude_only = 1;
    file_phase = [];
elseif isempty(file_phase)
    opt.flag_magnitude_only = 1;
else
    opt.flag_magnitude_only = 0;
end

if ~isfield(opt,'flag_discard_imag') || isempty(opt.flag_discard_imag)
    % When the low pass filter for phase is pretty narrow, keeping the
    % imaginary component becomes more important.
    opt.flag_discard_imag = 0;
end

if ~isfield(opt,'flag_force_no_mask') || isempty(opt.flag_force_no_mask)
    opt.flag_force_no_mask = 0;
end

% CB: Set phase normalization options
%   Defaults from original code were 3 and 'tukey', respectively
if ~isfield(opt,'phase_filter_width') || isempty(opt.phase_filter_width)
    opt.phase_filter_width = 4;
end
if ~isfield(opt,'phase_filter_width') || isempty(opt.phase_filter_width)
    opt.phase_filter_type = 'gauss';
end

%% check if bm4d files exist in the path

if exist('BM4D','file') ~= 2
    error(['BM4D files were not found.  ' ...
        '\n\nPlease download it at https://webpages.tuni.fi/foi/GCF-BM3D  (version 4.2.5 or %s)' ...
        '\nMake sure to extract the zip file into aux_functions.' ...
        '\n\nCheck Readme for further instructions.'],'newer');

end



%% *** prepare files, perform phase rotation if complex data ****
fprintf('\n **************************************** \n')


[I_noisy, I_noisy_imag] = dwi_prepareFiles(file_magn,file_phase,opt);

[filepath,~,~] = fileparts(file_magn) ;

mask_file = dir([filepath filesep '*mask_PSD.nii']); % mask for psd calculation. to avoid fake/forced backgrounds that are all zeroed out.

if ~isempty(mask_file) && ~opt.flag_force_no_mask
    I_mask_to_PSD = double(niftiread([mask_file(1).folder filesep mask_file(1).name]));
    I_mask_to_PSD = I_mask_to_PSD(:,:,:,1);
elseif opt.flag_force_no_mask
    sz = [size(I_noisy), 1, 1, 1];
    I_mask_to_PSD = ones(sz(1:3));
else
    % When empty, the mask will be automatically determined
    I_mask_to_PSD = [];
end


%% ** start the core processing ***

fprintf('Estimating noise Map and noise PSD ...\n')

bvals_file_path = file_magn;
bvals_file_path = erase(bvals_file_path,'.nii.gz');
bvals_file_path = erase(bvals_file_path,'.nii');
bvals_file_path = [bvals_file_path '.bval'];

bval_file = dir(bvals_file_path);

if ~isempty(bval_file)
    bvals = load(bvals_file_path);
    bvals_round = round(bvals./1000).*1000;
    uniq_bvals = unique(bvals_round);
else
    fprintf('No bval file found... using entire dataset for noise estimation ...\n')
    bvals_round = ones(1,size(I_noisy,4)); %fake bvals, to use all data for estimating noise
    uniq_bvals = 1;
end

save_phase = 0;
if opt.flag_magnitude_only == 0 %
  
    Nstage = 2;
    if opt.flag_discard_imag
        Nstage = 1;
    end
    fprintf('Stage 1 of %d. Estimating noise Map and noise PSD for real channel ...\n',Nstage)

    [Sigma_map,loc_PSD,I_mask_to_PSD] = estimate_noise_map_psd(I_noisy(:,:,:,bvals_round==max(uniq_bvals)),I_mask_to_PSD);
    [glob_PSD, PSD65] = upsamplePSD(Sigma_map, loc_PSD);

    [I_denoised] = BM4D_PC(I_noisy, Sigma_map, glob_PSD, PSD65, opt.profile_bm4d);

    if ~opt.flag_discard_imag
        fprintf('Stage 2 of %d. Estimating noise Map and noise PSD for imag channel ...\n',Nstage)
        % Note that here we re-use the I_mask_to_PSD estimated from real
        % data
        [Sigma_map_imag,loc_PSD_imag] = estimate_noise_map_psd(I_noisy_imag(:,:,:,bvals_round==max(uniq_bvals)),I_mask_to_PSD);
        [glob_PSD_imag, PSD65_imag] = upsamplePSD(Sigma_map_imag, loc_PSD_imag);
    
        [I_denoised_imag] = BM4D_PC(I_noisy_imag, Sigma_map_imag, glob_PSD_imag, PSD65_imag, opt.profile_bm4d);
        I_denoised_phase = atan2(I_denoised_imag, I_denoised);
        I_denoised = sqrt(I_denoised.^2 + I_denoised_imag.^2);
        clear I_denoised_imag

        save_phase = 1;
    else
        I_denoised = abs(I_denoised);
    end

else %magnitude data

    

    [Sigma_map,loc_PSD] = estimate_noise_map_psd_VST(I_noisy(:,:,:,bvals_round==max(uniq_bvals)),I_mask_to_PSD); 
    [glob_PSD, PSD65] = upsamplePSD(Sigma_map, loc_PSD);

    fprintf('Denoising Magnitude Image with VST\n')

    I_denoised = BM4D_PC_VST(I_noisy, Sigma_map, glob_PSD, PSD65, opt.profile_bm4d);
    I_denoised = max(0,I_denoised);


end


%% Save Denoised file 
dwi_savefile(file_magn,I_denoised,0)
if save_phase
    dwi_savefile(file_magn,I_denoised_phase,1)
end


fprintf('\nFinished denoising with BM4D-PC.\n')
fprintf('\n **************************************** \n')


end

function [glob_PSD, PSD65] = upsamplePSD(Sigma_map, loc_PSD)
    [X, Y, Z] = size(Sigma_map);

    if size(loc_PSD,3)>1
        glob_PSD = FFT_PSD_Upsampling(loc_PSD, [X Y Z]); % upsampling psd
        PSD65 = FFT_PSD_Upsampling(loc_PSD, [65 65 65]); % upsampling psd for bm4d parameters estimation
    else
        glob_PSD = repmat(FFT_PSD_Upsampling(loc_PSD, [X Y 1]),[1 1 Z]).*Z; % upsampling psd
        PSD65 = repmat(FFT_PSD_Upsampling(loc_PSD, [65 65 1]),[1 1 65]).*65; % upsampling psd for bm4d parameters estimation
    end
end
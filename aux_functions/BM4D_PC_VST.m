
function [I_PCA_denoised]= BM4D_PC_VST(I_Noisy, Estimated_NoiseMap, glob_PSD, PSD_65, profile)

if(~exist('profile', 'var'))
    profile = 'lc';
end

%****** Deviding noisy image by the Estimated map ****%
[X, Y, Z, W] = size(I_Noisy);

I_Noisy = I_Noisy./Estimated_NoiseMap;

%** Apply (forward) VST **

I_Noisy = riceVST(I_Noisy,1,'A');

% ****** PCA decomposition for denoising **********%

I_Noisy = reshape(I_Noisy,[X*Y*Z W]);

I_noisy_cov = I_Noisy'*I_Noisy;


[~, ~, v] = svd(I_noisy_cov, 'econ');


I_Noisy = I_Noisy*v;

I_Noisy = reshape(I_Noisy,[X Y Z W]);


%% ** adjust psd


glob_PSD = glob_PSD /  (norm(glob_PSD(:),1)/numel(glob_PSD)^2); % Just enforcing/guaranteeing it is normalized (var = 1)

variances = ones(1,W); % due to the normalization by Estimated map and VST, variance is one on every PC

%% bm4d profile setup
BM4D_profile= get_BM4D_profile(profile, PSD_65); %get profile settings

stage_arg = BM4DProfile.ALL_STAGES;


%% Denoising first Component and getting block-matching

fprintf('Denoising first Component and getting block-matching poisitions...\n')

I_PCA_denoised = I_Noisy;

i=1;

tic_first_pc = tic;

[I_PCA_denoised(:,:,:,i), match_arrs] = BM4D(I_Noisy(:,:,:,i), glob_PSD.*variances(i), BM4D_profile, stage_arg, {true, true});

elap_time_first_pc = toc(tic_first_pc);

fprintf('Elapsed time denoising first PC = %.2f (sec) ...\n',elap_time_first_pc)

%% Denoising remaining PCs
fprintf('Denoising remaining Components ...\n')

est_time_min = elap_time_first_pc * (W - 1) / 60;


fprintf('Remaining processing time = %.2f min ...\n\n', est_time_min);

finish_time = datetime('now') + minutes(est_time_min);
finish_time.Format = 'HH:mm';

fprintf('Estimated finish time: %s\n', string(finish_time));


for i=2:W
    [I_PCA_denoised(:,:,:,i)] = BM4D(I_Noisy(:,:,:,i), glob_PSD.*variances(i), BM4D_profile, stage_arg, match_arrs);
end

%%  'Inverse' PCA
I_PCA_denoised = reshape(I_PCA_denoised,[X*Y*Z W]);

I_PCA_denoised = I_PCA_denoised*v' ;

I_PCA_denoised = reshape(I_PCA_denoised,[X Y Z W]);


%  Inverse VST
I_PCA_denoised = riceVST_EUI(I_PCA_denoised,1,'A');


%******* Multiplying denoised image by the Estimated map ****%
I_PCA_denoised = I_PCA_denoised.*Estimated_NoiseMap;

%% mitigate zero values iside the brain due to very low snr and inverse vst 

I_PCA_denoised_aux = zeros(size(I_PCA_denoised));

template = fspecial3("average",3);

for i=1:size(I_PCA_denoised,4)

    I_PCA_denoised_aux(:,:,:,i) = imfilter(I_PCA_denoised(:,:,:,i),template,"symmetric");


end

I_PCA_denoised(I_PCA_denoised==0) = I_PCA_denoised_aux(I_PCA_denoised==0);


end

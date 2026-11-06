
function [Sigma_map, loc_PSD, I_mask_to_PSD]= estimate_noise_map_psd_VST(I_noisy_orig,I_mask_to_PSD,noisemap_filt_width)


% ****** PCA for noise estimation **********%

[X, Y, Z, W] = size(I_noisy_orig);

% estimate sigma map from last few PCs
num_of_PCs = 3;

% index of PC - backwards
idx = W - num_of_PCs +1;

[Sigma_map,loc_PSD,I_mask_to_PSD] = estimate_noise_map_psd(I_noisy_orig,I_mask_to_PSD);

if ~isempty(loc_PSD)
    [X, Y, Z] = size(Sigma_map);

    if size(loc_PSD,3)>1
        glob_PSD = FFT_PSD_Upsampling(loc_PSD, [X Y Z]); % upsampling psd
        PSD65 = FFT_PSD_Upsampling(loc_PSD, [65 65 65]); % upsampling psd for bm4d parameters estimation
    else
        glob_PSD = repmat(FFT_PSD_Upsampling(loc_PSD, [X Y 1]),[1 1 Z]).*Z; % upsampling psd
        PSD65 = repmat(FFT_PSD_Upsampling(loc_PSD, [65 65 1]),[1 1 65]).*65; % upsampling psd for bm4d parameters estimation
    end
end


I_denoised = BM4D_PC(I_noisy_orig, Sigma_map, glob_PSD, PSD65, 'mp', 1);

I_denoised = max(eps,I_denoised);

tempalte = fspecial3("average",7);

for i=1:size(I_denoised,4)
    denoised_aux = imfilter(I_denoised(:,:,:,i),tempalte,"symmetric");
    [~, Sigma_map_aux(:,:,:,i)] = ricePairInversion(denoised_aux,Sigma_map);
end

Sigma_map = mean(Sigma_map_aux,4);


%% Create mask for regions that may have been masked during recon, or if skull stripping was performed etc.
% We don't just check for 0, because iterative non-cartesian recons can
% have a tiny amount of signal leaking into the masked region (if receiver
% sensitivities are masked, for example).
if isempty(I_mask_to_PSD)
    % Should actually never enter this if statement, because
    % estimate_noise_map_psd above should find the mask.
    I_mask_to_PSD = estimateMask(Sigma_map);
end

% smooth noise map
Sigma_map = smooth_noise_map(Sigma_map,I_mask_to_PSD,noisemap_filt_width);

Sigma_map = max(eps,Sigma_map);


%% ***** estimate noise PSD in VST domain ******* %

I_noisy_orig = I_noisy_orig./Sigma_map;


I_noisy_orig = riceVST(I_noisy_orig,1,'B');


I_noisy = reshape(I_noisy_orig,[X*Y*Z W]);

[I_noisy, vals, ~] = svd(I_noisy, 'econ');

I_noisy = reshape(I_noisy,[X Y Z W]);

I_Noisy_to_est_PSD = I_noisy(:,:,:,idx:end);

singvals = diag(vals);
singvals = singvals(idx:end);


% ** only inside a masked region
I_Noisy_to_est_PSD = I_Noisy_to_est_PSD.*(repmat(I_mask_to_PSD,[1 1 1 size(I_Noisy_to_est_PSD,4)]));

[x1,y1,z1] = ind2sub(size(I_mask_to_PSD),find(I_mask_to_PSD));

I_Noisy_to_est_PSD = I_Noisy_to_est_PSD(unique(x1),unique(y1),unique(z1),:);

I_Noisy_to_est_PSD = I_Noisy_to_est_PSD(3:end-3,3:end-3,:,:); %avoid borders


% ** overlapping chunks of slices

chunk_size = 5;

overlap = 2;

% Calculate the step size for chunking (chunk_size - overlap)
step = chunk_size - overlap;

z_size = size(I_Noisy_to_est_PSD,3);


for i=1:num_of_PCs

    I_Noisy_to_est_PSD(:,:,:,i) = (I_Noisy_to_est_PSD(:,:,:,i).*singvals(i));


    counter = 0;

    for start_z=1:step:z_size

        % Define the end of the current chunk, making sure it does not exceed the image size
        end_z = min(start_z + chunk_size - 1, z_size);

        if abs(end_z - start_z+1)>3 %minimum of three slices. Avoid single slices on the end of the chunk

            counter = counter + 1;
            loc_PSD(:,:,counter,i)   = Est_rootPSD(I_Noisy_to_est_PSD(:,:,start_z:end_z,i)).^2;

        end
    end


end

loc_PSD = min(loc_PSD,[],3); % minimum projection intensity

loc_PSD = mean(loc_PSD,4);



end


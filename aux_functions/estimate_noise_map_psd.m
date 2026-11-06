
function [sigma_map, loc_PSD, I_mask_to_PSD]= estimate_noise_map_psd(I_noisy_orig,I_mask_to_PSD,noisemap_filt_width)

if nargin<3 || isempty(noisemap_filt_width)
    noisemap_filt_width = 4;
end


% ****** PCA for noise estimation **********%

[X, Y, Z, W] = size(I_noisy_orig);


% ****** PCA decomposition for denoising **********%

I_Noisy = reshape(I_noisy_orig,[X*Y*Z W]);

I_noisy_cov = I_Noisy'*I_Noisy;


[~, ~, v] = svd(I_noisy_cov, 'econ');


I_Noisy = I_Noisy*v;

I_Noisy = reshape(I_Noisy,[X Y Z W]);

%% estimate sigma map from last few PCs
num_of_PCs = 3;

% index of PC - backwards
idx = W - num_of_PCs +1;

I_noisy_to_est = I_Noisy(:,:,:,idx:end);


sigma_map = zeros([size(I_Noisy,1:3),num_of_PCs]);

for i=1:num_of_PCs

    template = fspecial3("average",5);

    E = imfilter(I_noisy_to_est(:,:,:,i) ,template,"symmetric"); %local mean (1st moment)
    E2 = imfilter(I_noisy_to_est(:,:,:,i) .^2,template,"symmetric");%local mean2 (2nd moment)


    sigma_map(:,:,:,i) = sqrt( max(0,E2 - E.^2) ); %var = E2 - E.^2

end

sigma_map = mean(sigma_map,4);

%% Create mask for regions that may have been masked during recon, or if skull stripping was performed etc.
% We don't just check for 0, because iterative non-cartesian recons can
% have a tiny amount of signal leaking into the masked region (if receiver
% sensitivities are masked, for example).

if isempty(I_mask_to_PSD)
    % Find the median with a cube with half the FOV to get typical sigma
    % values. Use central region to reduce how much masked regions affect the
    % median estimate. We also do some mild median filtering to lower the
    % chance of a random single voxel deep in the brain not reaching the
    % threshold for masking.
    sigma_map_a = medfilt3(sigma_map,3*[1 1 1]);
    sig_typ = sigma_map_a(round(end/4):round(end*3/4),round(end/4):round(end*3/4),round(end/4):round(end*3/4));
    sig_typ = median(sig_typ(:));
    
    % Set the threshold to 20% of median, since it would be very strange to
    % have the noise vary by so much.
    thresh_s = 0.2*sig_typ;
    I_mask_to_PSD =  sigma_map_a > thresh_s;
    clear sig_typ thresh_s

    frac_masked = 1 - sum(I_mask_to_PSD(:))/numel(I_mask_to_PSD);
    if frac_masked > 0.05
        fprintf('  Found that %.1f %% of voxels were likely masked. Not including them in PSD estimate.\n', 100*frac_masked )
    end
end

%% Smooth the noise map
% CB: median filter, implemented with medfilt3, has disadvantage of only cube shaped
% kernels. Median is advantageous over mean when there are large spikes,
% which does not seem to be likely for a noise map, so I will switch
% to gaussian smoothing with where 2*std = width supplied.
% sigma_map = medfilt3(sigma_map,med_filt_width*[1 1 1]);
alph = (noisemap_filt_width-1)/noisemap_filt_width;
gw = gausswin(2*noisemap_filt_width+1,2*alph); 
kern = gw(:).*reshape(gw,[1 2*noisemap_filt_width+1]).*reshape(gw,[1 1 2*noisemap_filt_width+1]);
kern = kern/sum(kern(:));

% Use a weighted convolution to suppress contributions of zeros from
% outside the mask
w = convn(I_mask_to_PSD,kern,'same')+0.0001; % weights for weighted convolution. Add a small value to prevent div by 0
sigma_map = convn(sigma_map.*w,kern,'same')./convn(w,kern,'same');

sigma_map(sigma_map<=0)=eps;


%% *********** estimate noise PSD, complex-valued data ***********

I_Noisy_to_est_PSD = I_noisy_to_est./sigma_map; %normalize


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

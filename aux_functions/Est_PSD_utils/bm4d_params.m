function [params] = bm4d_params(psd)

params.lambda = [];
params.lambda_re = [];
params.mu = [];
params.mu_re = [];


[params.lambda, params.mu2, params.lambda_re, params.mu2_re] = process_psd(psd);



end


function [lambda, wielambdasq, lambda2, wielambdasq2] = estimate_psd_parameters(PSD65_full)
% Estimate parameters based on the sigma_psd

% Get the optimal parameters and matching features for a bunch of PSDs
load('param_matching_data.mat', 'features', 'maxes');

sz = 65;
indices_to_take = [1:2:10 12:5:32];

lambda = [];
wielambdasq = [];
lambda2 = [];
wielambdasq2 = [];

% Get separate parameters for each sigma_psd provided
for psd_num = 1:size(PSD65_full, 3)
    PSD65 = PSD65_full(:, :, psd_num);
    PSD65 = fftshift(PSD65);
    
    % Get features for this sigma_psd
    pcaxa = get_features(PSD65, sz, indices_to_take);

    % Calculate distances to other PSDs
    mm = mean(features, 2);
    centered_features = features - mm;
    corr_matx = centered_features * centered_features';
    corr_matx = corr_matx / 500;
    
    centered_pcax = pcaxa' - mm;
    
    [u, s, ~] = svd(corr_matx);
    centered_features = u * centered_features;
    centered_pcax = u * centered_pcax;
    centered_features = centered_features .* sqrt(diag(s));
    centered_pcax = centered_pcax .* sqrt(diag(s));

    diff_pcax = sqrt(sum(abs(centered_features - centered_pcax).^2, 1));

    % Take only smallest->best x %
    [~, dff_I] = sort(diff_pcax);

    %  Take 20 most similar PSDs into consideration
    count = 20;
    diff_indices = dff_I(1:count);

    % Invert, smaller -> bigger weight
    diff_inv = 1 ./ (diff_pcax + eps);
    diff_inv = diff_inv(diff_indices) ./ sum(diff_inv(diff_indices));

    % Weight
    param_idxs = sum(diff_inv .* maxes(diff_indices, :)', 2);

    lambdas = 2.5:0.1:5;
    wielambdasqs = 0.2:0.2:6;
    
    % Get parameters from indices - 
    % Interpolate lambdas and mu^2s from the list
    for ix = [1, 3]
        param_idx = max(1, param_idxs(ix));
        param_idx2 = max(1, param_idxs(ix+1));

        l1 = lambdas(floor(param_idx));
        l2 = lambdas(min(ceil(param_idx), numel(lambdas)));

        w1 = wielambdasqs(floor(param_idx2));
        w2 = wielambdasqs(min(ceil(param_idx2), numel(wielambdasqs)));

        param_smooth = param_idx - floor(param_idx);
        param_smooth2 = param_idx2 - floor(param_idx2);

        if ix == 1
            lambda = [lambda, l2 * param_smooth + l1 * (1 - param_smooth)];
            wielambdasq = [wielambdasq, w2 * param_smooth2 + w1 * (1 - param_smooth2)];
        elseif ix == 3
            lambda2 = [lambda2, l2 * param_smooth + l1 * (1 - param_smooth)];
            wielambdasq2 = [wielambdasq2, w2 * param_smooth2 + w1 * (1 - param_smooth2)];
        end

    end
end
end

% Calculate features for a sigma_psd from integrals
function f = get_features(sigma_psd, sz, indices_to_take)

    [I_rot, I_rot2] = pcax(sigma_psd);
    f1 = zeros(1, numel(indices_to_take));
    f2 = f1;
    
    % Extract features for a sigma_psd
    for ii = 1:numel(indices_to_take)
        rang = indices_to_take(ii);
        if ii > 1
            rang = indices_to_take(ii-1)+1:rang;
        end
        f1(ii) = sum(I_rot(ceil(sz/2) + rang - 1)) / numel(rang);
        f2(ii) = sum(I_rot2(ceil(sz/2) + rang - 1)) / numel(rang);
    end
    
    f = [f1 f2];
end

% Calculate integrals along principal axes of the sigma_psd
function [I_rot, I_rot2] = pcax(sigma_psd)

N=size(sigma_psd,1);
[G2, G1]=meshgrid(1:N,1:N);

trapz2D=@(G2,G1,sigma_psd) trapz(G1(:,1),trapz2(G2,sigma_psd,2),1);

Pn=sigma_psd/trapz2D(G2,G1,sigma_psd);

m2=trapz2D(G2,G1,Pn.*G2);
m1=trapz2D(G2,G1,Pn.*G1);
C=zeros(2);
O1=[2 1 1 0];
O2=[0 1 1 2];
for jj=[1 2 4]
    C(jj)=squeeze(trapz2D(G2,G1,Pn.*(G2-m2).^O1(jj).*(G1-m1).^O2(jj)));
end
C(3)=C(2);


[U, ~, ~]=svd(C);

N3 = 3 * N;
[G13N, G23N]=ndgrid((1:N3)-(N3+1)/2, (1:N3)-(N3+1)/2);

% Rotate PSDs and calculate integrals along the rotated PSDs
theta = angle(U(1, 1) + 1i * U(1, 2));
G2rot=G23N(N+1:2*N,N+1:2*N)*cos(theta)-G13N(N+1:2*N,N+1:2*N)*sin(theta);
G1rot=G13N(N+1:2*N,N+1:2*N)*cos(theta)+G23N(N+1:2*N,N+1:2*N)*sin(theta);
P_rot_handle = griddedInterpolant(G13N,G23N,repmat(sigma_psd,[3,3]),'linear','nearest');
P_rot = P_rot_handle(G1rot,G2rot);
I_rot = trapz2(G1, P_rot, 1);

theta2 = angle(U(2, 1) + 1i * U(2, 2));
G2rot=G23N(N+1:2*N,N+1:2*N)*cos(theta2)-G13N(N+1:2*N,N+1:2*N)*sin(theta2);
G1rot=G13N(N+1:2*N,N+1:2*N)*cos(theta2)+G23N(N+1:2*N,N+1:2*N)*sin(theta2);
P_rot2 = P_rot_handle(G1rot,G2rot);
I_rot2 = trapz2(G1, P_rot2, 1);
end

function I = trapz2(X,Y,dimm)

if dimm==2
I=sum((Y(:,2:end)+Y(:,1:end-1))/2.*(X(:,2:end)-X(:,1:end-1)),2);
else
I=sum((Y(2:end,:)+Y(1:end-1,:))/2.*(X(2:end,:)-X(1:end-1,:)),1);
end

end

% Process sigma_psd, get parameters if needed
function [lambda_thr, mu2, lambda_thr_re, mu2_re] = process_psd(sigma_psd)

auto_params = true;


% Calculate the correlation kernel from the sigma_psd in order
% to resize the sigma_psd. (skip if we are not resizing the sigma_psd)
if (auto_params)
    temp_kernel = real(fftshift(ifftn(sqrt(sigma_psd / numel(sigma_psd)))));
end

if auto_params
    % Estimate parameters based on the sigma_psd
    
    minus_size = ceil((size(sigma_psd) - 65) / 2);
    temp_kernel65 = temp_kernel;
    if minus_size(1) > 0
        temp_kernel65 = temp_kernel65(1+minus_size(1):minus_size(1)+65, :, :);
    end
    if minus_size(2) > 0
        temp_kernel65 = temp_kernel65(:, 1+minus_size(2):minus_size(2)+65, :);
    end
    if numel(minus_size) > 3 && minus_size(3) > 0
        temp_kernel65 = temp_kernel65(:, :, 1+minus_size(3):minus_size(3)+65);
    end
    temp_kernel65 = temp_kernel65 / norm(temp_kernel65(:));
    if size(temp_kernel65, 3) > 1
        PSD65 = abs(fftn(temp_kernel65, [65, 65, 65])).^2 * 65 * 65;
        % Reduce PSD to two dimensions to use BM3D parameter estimation
        PSD65 = project_psd_for_param_est(fftshift(PSD65));
    else
        PSD65 = abs(fftn(temp_kernel65, [65, 65])).^2 * 65 * 65;
    end
    PSD65 = PSD65 / mean(PSD65(:)) * 65 * 65;

    % Parameter estimation
    [lambda_thr, mu2, lambda_thr_re, mu2_re] = estimate_psd_parameters(PSD65);
end

end


function P_rot_f = project_psd_for_param_est(sigma_psd)
% Return the 2-D projection of maximum variance squeezed along a principal
% axis.
% sigma_psd should be a NxNxN matrix.

N=size(sigma_psd,1);
[G2, G1, G3]=meshgrid(1:N,1:N,1:N);

trapz2D=@(G3,G2,G1,sigma_psd) trapz3(G2(:,:,1),trapz3(G3,sigma_psd,3),2);
trapz3D=@(G3,G2,G1,sigma_psd) trapz(G1(:,1,1), trapz2D(G3, G2, G1, sigma_psd)); 

Pn=sigma_psd/trapz3D(G3,G2,G1,sigma_psd);

m2=trapz3D(G3,G2,G1,Pn.*G2);
m1=trapz3D(G3,G2,G1,Pn.*G1);
m3=trapz3D(G3,G2,G1,Pn.*G3);

C=zeros(3);
O1=[2 1 1 1 0 0 1 0 0];
O2=[0 1 0 1 2 1 0 1 0];
O3=[0 0 1 0 0 1 1 1 2];
for jj=[1, 2, 3, 5, 6, 9]
    C(jj)=squeeze(trapz3D(G3,G2,G1,Pn.*(G2-m2).^O2(jj).*(G1-m1).^O1(jj).*(G3-m3).^O3(jj)));
end
C(4) = C(2);
C(7) = C(3);
C(8) = C(6);

% Select highest variance
[U, ~, ~]=svd(C);
P_rot_f = 0;
P_rot_var_f = 0;

for ax_ix=1:3
    v = (U(:, ax_ix));
    P_rot = ifftshift(rotateTo(sigma_psd, v));
    P_rot_proj = squeeze(P_rot(1, :, :));
    P_var = sum(P_rot_proj(:));
    if P_var > P_rot_var_f
        P_rot_f = P_rot_proj;
        P_rot_var_f = P_var;
    end
end
end

function P_rot = rotateTo(sigma_psd, v)
% Rotate sigma_psd to be oriented by v (3x1)
% sigma_psd should be a NxNxN matrix;
% v should be unit length.

r = [1, 0, 0];
if sum(abs(v(:)-r(:))) < 1e-3 || sum(abs(v(:)+r(:))) < 1e-3
    P_rot = sigma_psd;
    return;
end
rotAx = cross(v, r);
rotAx = rotAx / norm(rotAx(:));
rotAngle = acos(min(max(dot(v, r), -1), 1));

rotAxCr = [0, -rotAx(3), rotAx(2); rotAx(3), 0, -rotAx(1); -rotAx(2), rotAx(1), 0];
rotationMatrix = cos(rotAngle) * eye(3) + sin(rotAngle) * rotAxCr + ...
    (1 - cos(rotAngle)) * (rotAx' * rotAx);

N = size(sigma_psd, 1);
N3=3*N;
[G13N, G23N, G33N]=ndgrid((1:N3)-(N3+1)/2, (1:N3)-(N3+1)/2, (1:N3)-(N3+1)/2);
gns = cat(4, G13N, G23N, G33N);
Grot=zeros(N, N, N, 3);
for x=1:N
    for y=1:N
        for z=1:N
        rot = rotationMatrix * squeeze(gns(N+z,N+y,N+x, :));
        Grot(z, y, x, :) = rot;
        end
    end
end
P_rot_handle = griddedInterpolant(G13N,G23N,G33N,repmat(sigma_psd,[3,3,3]),'linear','nearest');
P_rot = P_rot_handle(Grot(:, :, :, 1), Grot(:, :, :, 2), Grot(:, :, :, 3));

end

function I = trapz3(X,Y,dimm)

if dimm==3
    I=sum((Y(:,:,2:end)+Y(:,:,1:end-1))/2.*(X(:,:,2:end)-X(:,:,1:end-1)),3);
elseif dimm==2
    I=sum((Y(:,2:end,:)+Y(:,1:end-1,:))/2.*(X(:,2:end,:)-X(:,1:end-1,:)),2);
else
    I=sum((Y(2:end,:,:)+Y(1:end-1,:,:))/2.*(X(2:end,:,:)-X(1:end-1,:,:)),1);
end

end
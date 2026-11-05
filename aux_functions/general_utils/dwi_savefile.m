
function dwi_savefile(file_magn,I_denoised,is_phase)

fprintf('Saving file ... \n');

save_file_prefix = 'BM4D_PC_';

file_dir = dir(file_magn);
file_folder = file_dir.folder;


fname_tosave = [save_file_prefix file_dir.name];
fname_tosave = erase(fname_tosave,'.nii.gz');
fname_tosave = erase(fname_tosave,'.nii');

if is_phase
    fname_tosave = [fname_tosave, '_phase'];
end

nii_info = niftiinfo([file_folder filesep file_dir.name]);

% Always save as at least single, but save as double if input data was
% double
if nii_info.raw.bitpix > 32
    I_denoised = cast(I_denoised,'double');
    nii_info.raw.datatype = 64;
    nii_info.raw.bitpix = 64;
    nii_info.Datatype = 'double';
    nii_info.BitsPerPixel = 64;
else
    I_denoised = cast(I_denoised,'single');
    nii_info.raw.datatype = 16;
    nii_info.raw.bitpix = 32;
    nii_info.Datatype = 'single';
    nii_info.BitsPerPixel = 32;
end

if contains(file_dir.name,'nii.gz')
    doCompress = true;
else
    doCompress = false;
end
niftiwrite(I_denoised,[file_folder filesep fname_tosave], nii_info, 'Compressed', doCompress);
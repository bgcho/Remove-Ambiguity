% run_remove_ambiguity.m
% Wrapper to remove the ambiguity with different options

% Load the images, Coordinates, and Bathymetry
load('../raw_data/edited_images.mat')
load('../raw_data/Lofoten_bathy.mat')
load('../raw_data/jet_bg.mat')

% Mask the region near source
near_src_range = 6300;
% Threshold the image
th = 138;
stride_all = [1 2 3];

for ii = 1:length(stride_all)
	stride = stride_all(ii);
	run('remove_ambiguity.m')
	clear all
end


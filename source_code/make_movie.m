function make_movie(im_all, timestamp, grid_x, grid_y, bathy_map, bathygrid_x, bathygrid_y, center_x, center_y, avi_name, frame_rate)

load('jet_bg.mat')
vidObj = VideoWriter(avi_name);
vidObj.FrameRate = frame_rate;
open(vidObj);


no_images = size(im_all,3);
figure; hold on
imagesc((grid_x-center_x)/1000, (grid_y-center_y)/1000, im_all(:,:,1))
[c_bathy, h_bathy] = contour((bathygrid_x-center_x)/1000, (bathygrid_y-center_y)/1000, bathy_map, [-300:50:0],'k--');
axis xy equal tight
colorbar
caxis([80 180])
title(datestr(timestamp(1)))
xlim([min((grid_x-center_x)/1000) max((grid_x-center_x)/1000)])
ylim([min((grid_y-center_y)/1000) max((grid_y-center_y)/1000)])
xlabel('Easting (km)')
ylabel('Northing (km)')
colormap(cmap)
set(gca,'nextplot','replacechildren');



for ii = 1:no_images*frame_rate
	frame_no = ceil(ii/frame_rate);
	imagesc((grid_x-center_x)/1000, (grid_y-center_y)/1000, im_all(:,:,frame_no)); hold on
	[c_bathy, h_bathy] = contour((bathygrid_x-center_x)/1000, (bathygrid_y-center_y)/1000, bathy_map, [-300:50:0],'k--');
	axis xy equal tight
	caxis([80 180])
	title(datestr(timestamp(frame_no)))
	xlim([min((grid_x-center_x)/1000) max((grid_x-center_x)/1000)])
	ylim([min((grid_y-center_y)/1000) max((grid_y-center_y)/1000)])

	% Write each frame to the file.
	currFrame = getframe;
	writeVideo(vidObj,currFrame);
end

% Close the file.
close(vidObj);


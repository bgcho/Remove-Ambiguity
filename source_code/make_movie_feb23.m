load('edited_images.mat')
load('Lofoten_bathy.mat')
load('jet_bg.mat')

no_images = size(im_all,3);
for ii = 1:no_images
	figure; hold on
	imagesc((grid_x-center_x)/1000, (grid_y-center_y)/1000, im_all(:,:,ii))
	[c_bathy, h_bathy] = contour((bathygrid_x-center_x)/1000, (bathygrid_y-center_y)/1000, bathy_map, [-300:50:0],'k--');
	axis xy equal tight
	colorbar
	caxis([80 180])
	title(datestr(timestamp(ii)))
	xlim([min((grid_x-center_x)/1000) max((grid_x-center_x)/1000)])
	ylim([min((grid_y-center_y)/1000) max((grid_y-center_y)/1000)])
	xlabel('Easting (km)')
	ylabel('Northing (km)')
	colormap(cmap)
	saveas(gcf, ['../result_images/' ping_id(ii,:)],'tif')
	saveas(gcf, ['../result_images/' ping_id(ii,:)],'epsc')
	saveas(gcf, ['../result_images/' ping_id(ii,:)],'png')
	saveas(gcf, ['../result_images/' ping_id(ii,:)],'fig')
	close all
end

run('make_movie_fromtif.m');




%save_file = '../result_images/Rost_feb23.avi';
%frame_rate = 5;

%make_movie(im_all, timestamp, grid_x, grid_y, bathy_map, bathygrid_x, bathygrid_y, center_x, center_y, save_file, frame_rate)


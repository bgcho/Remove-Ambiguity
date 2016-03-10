
frame_per_sec = 15;
    
tif_files = dir('../result_images/afambcorr_LK*.png');

M = moviein(length(tif_files));

filename = '../result_images/Rost_feb23_afambcorr_LK.avi';
aviObj = VideoWriter(filename);
aviObj.FrameRate = frame_per_sec;
open(aviObj);


for ii = 1:length(tif_files)*frame_per_sec
	file_count = ceil(ii/frame_per_sec);
	tif_filename = ['../result_images/' tif_files(file_count).name];

	imshow(tif_filename,'InitialMagnification','fit')
	axis equal tight
	M(ii).cdata = getframe(gca);
	writeVideo(aviObj, M(ii).cdata);
end
close(aviObj)



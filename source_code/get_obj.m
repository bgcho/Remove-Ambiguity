% Get the shoal object

function [contour_obj, mask_obj, cntr_obj, real_ghost_pair] = get_obj(grid_x, grid_y, im_all, src_utm_x, src_utm_y, heading, th, near_src_range)

tran_vec = [cosd(90-heading) ; sind(90-heading)];
norm_vec = [tran_vec(2,:) ; -tran_vec(1,:)];
for ii = 1:size(im_all,3)
	im_all(:,:,ii) = imfilter(im_all(:,:,ii),fspecial('gaussian',28,4),'same','replicate');
	% im_all(:,:,ii) = imfilter(im_all(:,:,ii),fspecial('disk',8),'same','replicate');
end



% mask_obj = {};
for ii = 1:size(im_all,3)
	th_mask(:,:,ii) = im_all(:,:,ii)>th;
	
	range_x = grid_x-src_utm_x(ii);
	range_y = grid_y-src_utm_y(ii);
	[Rx Ry] = meshgrid(range_x, range_y);
	Range = sqrt(Rx.^2+Ry.^2);
	far_mask(:,:,ii) = Range > near_src_range;
	object_mask(:,:,ii) = th_mask(:,:,ii) & far_mask(:,:,ii);
	% Get the contours and mask of each object
	[c_obj{ii}, ~] = contour(grid_x, grid_y, object_mask(:,:,ii), [1 1], 'm', 'linewidth',3);
	contour_indx{ii} = 1;
	while contour_indx{ii}(end)<size(c_obj{ii}, 2)
		contour_indx{ii} = [contour_indx{ii} contour_indx{ii}(end)+c_obj{ii}(2,contour_indx{ii}(end))+1];
	end
	contour_indx{ii} = contour_indx{ii}(1:end-1);
	contour_length{ii} = c_obj{ii}(2, contour_indx{ii});
	tmp_no_object = length(contour_indx{ii});
	jj = 1;	kk = 1;
	for jj = 1:tmp_no_object
		if jj < tmp_no_object
		tmp_contour_obj = c_obj{ii}(:,contour_indx{ii}(jj)+1:contour_indx{ii}(jj+1)-1);
		else
		tmp_contour_obj = c_obj{ii}(:,contour_indx{ii}(jj)+1:end);
		end
		if sum(tmp_contour_obj(:,end) ~= tmp_contour_obj(:,1))>0
			tmp_contour_obj = [tmp_contour_obj tmp_contour_obj(:,1)];
		end

		BW = roipoly(grid_x, grid_y, im_all(:,:,ii), tmp_contour_obj(1,:), tmp_contour_obj(2,:));
		if numel(find(BW>0))>100
			contour_obj{ii,kk} = tmp_contour_obj;
			mask_obj{ii}(:,:,kk) = BW;
			cntr_obj{ii}(1,kk) = sum(sum(BW.*(Rx+src_utm_x(ii))))/sum(sum(BW));
			cntr_obj{ii}(2,kk) = sum(sum(BW.*(Ry+src_utm_y(ii))))/sum(sum(BW));
			kk = kk+1;
		end
		
		% [lat,lon] = utm2deg(contour_shoal{ii,jj}(1,:).', contour_shoal{ii,jj}(2,:).', repmat(zone,[size(contour_shoal{ii,jj},2),1]));
		% contour_shoal_latlon{ii,jj} = [lon.' ; lat.'];
	end

	% Calculate the axis-symmetric position of each object center
	no_object = size(mask_obj{ii}, 3);
	src_utm = repmat([src_utm_x(ii) ; src_utm_y(ii)], [1,no_object]);
	ax_sym_cntr_obj{ii} = src_utm + tran_vec(:,ii)*(tran_vec(:,ii).'*(cntr_obj{ii}-src_utm)) - norm_vec(:,ii)*(norm_vec(:,ii).'*(cntr_obj{ii}-src_utm));

	% Pair the real and ambiguous objects
	for ll = 1:no_object
		[~, amb_indx] = min((ax_sym_cntr_obj{ii}(1,:)-cntr_obj{ii}(1,ll)).^2 + (ax_sym_cntr_obj{ii}(2,:)-cntr_obj{ii}(2,ll)).^2);
		real_ghost_pair{ii}(:,ll) = [ll ; amb_indx];		
	end

end

close all




% Calculate the SIFT flow of each image and get the average SIFT flow within each object

% Compare between the pairs

% Remove the ambiguous data





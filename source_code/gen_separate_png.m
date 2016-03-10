load('9_amb_corr.mat')
flow_type = {'LK'};
for yy = 1:length(flow_type)
	type = flow_type{yy};

	switch lower(type)
	case 'lk'
		flow_mag = LKflow_mag;
	case 'hs'
		flow_mag = HSflow_mag;
	case 'sift'
		flow_mag = SIFTflow_mag;
	end

	clear real_obj ghost_obj

	for ii = 1:no_image
		real_obj{ii} = [];
		ghost_obj{ii} = [];
		no_object = size(real_ghost_pair{ii},2);
		for jj = 1:no_object
			flow_obj1 = flow_mag(:,:,ii).*mask_obj{ii}(:,:,real_ghost_pair{ii}(1,jj));
			avg_flow_obj1 = sum(sum(flow_obj1))/numel(find(flow_obj1>0));
			std_flow_obj1 = sqrt(sum(sum(flow_obj1.^2))/numel(find(flow_obj1>0)) - avg_flow_obj1^2);
			flow_obj2 = flow_mag(:,:,ii).*mask_obj{ii}(:,:,real_ghost_pair{ii}(2,jj));
			avg_flow_obj2 = sum(sum(flow_obj2))/numel(find(flow_obj2>0));
			std_flow_obj2 = sqrt(sum(sum(flow_obj2.^2))/numel(find(flow_obj2>0)) - avg_flow_obj2^2);
			if avg_flow_obj2-avg_flow_obj1>max(std_flow_obj1, std_flow_obj2)
				real_obj{ii} = [real_obj{ii} real_ghost_pair{ii}(1,jj)];
				ghost_obj{ii} = [ghost_obj{ii} real_ghost_pair{ii}(2,jj)];
			elseif abs(avg_flow_obj1-avg_flow_obj2)<=max(std_flow_obj1, std_flow_obj2)
				real_obj{ii} = [real_obj{ii} real_ghost_pair{ii}(1,jj)];
			else
				real_obj{ii} = [real_obj{ii} real_ghost_pair{ii}(2,jj)];
				ghost_obj{ii} = [ghost_obj{ii} real_ghost_pair{ii}(1,jj)];
			end
		end
		real_obj{ii} = unique(real_obj{ii});
		ghost_obj{ii} = unique(ghost_obj{ii});
	end

	% Remove the data in the ambiguous object 
	mask_ambig = false(length(grid_y), length(grid_x), no_image);
	for ii = 1:no_image
		for jj = 1:length(ghost_obj{ii})
			mask_ambig(:,:,ii) = mask_ambig(:,:,ii) | mask_obj{ii}(:,:,ghost_obj{ii}(jj));		
		end	
		% neighbor = imfilter(double(mask_ambig(:,:,ii)), fspecial('gaussian',32,8), 'same','replicate')>0.01 & ~mask_ambig(:,:,ii);
		mask_ambig(:,:,ii) = imfilter(double(mask_ambig(:,:,ii)), fspecial('disk',30), 'same','replicate')>2e-04;
		neighbor = imfilter(double(mask_ambig(:,:,ii)), fspecial('disk',30), 'same','replicate')>4e-04 & ~mask_ambig(:,:,ii);
		mask_neighbor(:,:,ii) = neighbor;
		cur_image = im_all(:,:,ii);
		% Get the level of the neighbor region of the ambiguity
		level_neighbor{ii} = cur_image(find(neighbor>0));
		[N(:,ii), X(:,ii)] = hist(level_neighbor{ii}, 100);
	end

	im_ambig_rm = im_all.*(~mask_ambig);

	% Fill in the ambiguity region with data from the neighbor region
	% Get the pdf of the neighbor
	glob_level = linspace(min(min(X)), max(max(X)),100).';
	d_level = glob_level(2)-glob_level(1);

	for ii = 1:no_image
		glob_count(:,ii) = interp1(X(:,ii), N(:,ii), glob_level,'nearest','extrap');
		norm_factor = sum(glob_count(:,ii))*d_level;
		glob_pdf(:,ii) = glob_count(:,ii)/norm_factor;
	end

	% Fill pseudo data in the ambiguity region
	for ii = 1:no_image
		lin_indx = find(mask_ambig(:,:,ii)>0);
		if ~isempty(lin_indx)
		rand_num = gen_rand_pdf(glob_level, glob_pdf(:,ii), numel(lin_indx));
		cur_image = im_all(:,:,ii);
		cur_image(lin_indx) = rand_num;
		end
		im_pseudo(:,:,ii) = cur_image;
	end


	for ii = 1:no_image
		figure
		imagesc((grid_x-center_x)/1000, (grid_y-center_y)/1000, im_pseudo(:,:,ii)); hold on
		[c_bathy, h_bathy] = contour((bathygrid_x-center_x)/1000, (bathygrid_y-center_y)/1000, bathy_map, [-300:50:0],'k--');
		axis xy equal tight
		caxis([80 170])
		set(gca,'fontsize',12, 'fontweight','normal')
		title(datestr(timestamp(ii)))
		xlim([min((grid_x-center_x)/1000) max((grid_x-center_x)/1000)])
		ylim([min((grid_y-center_y)/1000) max((grid_y-center_y)/1000)])
		xlabel('Easting (km)'); ylabel('Northing (km)'); cbh = colorbar; set(cbh,'fontsize',12,'fontweight','normal'); ylabel(cbh, 'Intensity')
		colormap(cmap)
		saveas(gcf,['../result_images/afambcorr_LK_' num2str(ii)],'png')
		close
	end
	
end

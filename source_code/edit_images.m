% Cut the images at Lofoten to appropriate size

glob_xmax = 998952.8583969474;
glob_xmin = 962712.3320811579;
glob_xcnt = (glob_xmax+glob_xmin)/2;
glob_ymax = 7506121.700017895;
glob_ymin = 7468437.489491579;
glob_ycnt = (glob_ymax+glob_ymin)/2;
%grid_inc = 30;
no_glob_xgrid = round((glob_xmax - glob_xmin)/30);
no_glob_ygrid = round((glob_ymax - glob_ymin)/30);
glob_x = linspace(glob_xmin, glob_xmax, no_glob_xgrid);
glob_y = linspace(glob_ymin, glob_ymax, no_glob_ygrid);

% gfilter = fspecial('gaussian', 16, 4);

eval(['!ls ' '../raw_images/*.arr > rawimage_list'])
eval(['!ls ' '../raw_images/*.m > grid_list'])
eval(['!ls ' '../raw_images/bimap* > bimap_list'])

fid_arr = fopen('rawimage_list','r');
fid_grid = fopen('grid_list','r');
fid_bimap = fopen('bimap_list','r');
kk = 1;
while ~feof(fid_arr)
	arr_name = fgetl(fid_arr);
	grid_name = fgetl(fid_grid);
	im = read_arr(arr_name);
	im = flipud(im);
	im = (im-min(min(im)))*255/(max(max(im))-min(min(im)));
	run(grid_name);
	grid_x = [grid_xmin:grid_inc:grid_xmax];
	grid_y = [grid_ymin:grid_inc:grid_ymax];
	[~, xcnt_indx] = min(abs(glob_xcnt-grid_x));
	[~, ycnt_indx] = min(abs(glob_ycnt-grid_y));
	indx_x = [xcnt_indx-floor((no_glob_xgrid-1)/2):xcnt_indx+ceil((no_glob_xgrid-1)/2)];
	indx_y = [ycnt_indx-floor((no_glob_ygrid-1)/2):ycnt_indx+ceil((no_glob_ygrid-1)/2)];
	im = im(indx_y,indx_x);
	im_all(:,:,kk) = im;
	% im_all(:,:,kk) = 10*log10(conv2(10.^(im/10), gfilter,'same'));
	
	% Source locations of the pings
	bimap_name = fgetl(fid_bimap);
	fid_bimap_each = fopen(bimap_name,'r');
	mm = 1;
	while ~feof(fid_bimap_each)
		entry = fgetl(fid_bimap_each);
		delim = regexp(entry, '=');
		bimap_info{mm} = entry(delim+2:end-1);
		mm = mm+1;
	end
	src_lon(kk) = str2num(bimap_info{2});
	src_lat(kk) = str2num(bimap_info{3});
	src_utm_x(kk) = str2num(bimap_info{4});
	src_utm_y(kk) = str2num(bimap_info{5});
	heading(kk) = str2num(bimap_info{14});
	
	% Time stamp of the pings
	delim = regexp(arr_name,'/');
	ping_id = arr_name(delim(end)+1:end-4);
	ping_hour = str2num(ping_id(end-5:end-4));
	ping_min = str2num(ping_id(end-3:end-2));
	ping_sec = str2num(ping_id(end-1:end));
	jdcount = str2num(ping_id(end-9:end-7));
	ping_date = JDN2Greg(2014, jdcount);
	timestamp(kk) = datenum([ping_date ping_hour ping_min ping_sec]);
	ping_id_all(kk,:) = ping_id;
	
	kk = kk+1;
end
center_x = glob_xcnt;
center_y = glob_ycnt;
grid_x = glob_x;
grid_y = glob_y;
ping_id = ping_id_all;

save('../raw_images/edited_images.mat','ping_id','grid_x','grid_y','center_x', 'center_y','im_all', 'timestamp','src_lon','src_lat','src_utm_x','src_utm_y','heading');


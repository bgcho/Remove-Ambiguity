
%im1 = im2double(rgb2gray(imread('/home/bgcho/CourseWork/6_869/Pset3AB/Pset3AB/ps3_motion/car1.jpg')));
%im2 = im2double(rgb2gray(imread('/home/bgcho/CourseWork/6_869/Pset3AB/Pset3AB/ps3_motion/car2.jpg')));
glob_xmax = 998922.8583969474;
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


Sx = [-1 1 ; -1 1];
Sy = [-1 -1 ; 1 1];
gfilter = fspecial('gaussian', 16, 4);


im1 = read_arr('/home/bgcho/CourseWork/6_869/Project/raw_images/fora2014jd054t143319.arr');
im1 = flipud(im1);
%im1 = (im1-min(min(im1)))*255/(max(max(im1))-min(min(im1)));
run('/home/bgcho/CourseWork/6_869/Project/raw_images/fora2014jd054t143319.m')
grid_x1 = [grid_xmin:grid_inc:grid_xmax];
grid_y1 = [grid_ymin:grid_inc:grid_ymax];
[~, xcnt_indx1] = min(abs(glob_xcnt-grid_x1));
[~, ycnt_indx1] = min(abs(glob_ycnt-grid_y1));
indx_x1 = [xcnt_indx1-floor((no_glob_xgrid-1)/2):xcnt_indx1+ceil((no_glob_xgrid-1)/2)];
indx_y1 = [ycnt_indx1-floor((no_glob_ygrid-1)/2):ycnt_indx1+ceil((no_glob_ygrid-1)/2)];
im1 = im1(indx_x1, indx_y1);
im1 = 10*log10(conv2(10.^(im1/10), gfilter,'same'));

im1_x = conv2(10.^(im1/10), Sx, 'same');
im1_y = conv2(10.^(im1/10), Sy, 'same');
im1_grad = 10*log10(sqrt(im1_x.^2+im1_y.^2));

im2 = read_arr('/home/bgcho/CourseWork/6_869/Project/raw_images/fora2014jd054t143459.arr');
im2 = flipud(im2);
%im2 = (im2-min(min(im2)))*255/(max(max(im2))-min(min(im2)));
run('/home/bgcho/CourseWork/6_869/Project/raw_images/fora2014jd054t143459.m')
grid_x2 = [grid_xmin:grid_inc:grid_xmax];
grid_y2 = [grid_ymin:grid_inc:grid_ymax];
[~, xcnt_indx2] = min(abs(glob_xcnt-grid_x2));
[~, ycnt_indx2] = min(abs(glob_ycnt-grid_y2));
indx_x2 = [xcnt_indx2-floor((no_glob_xgrid-1)/2):xcnt_indx2+ceil((no_glob_xgrid-1)/2)];
indx_y2 = [ycnt_indx2-floor((no_glob_ygrid-1)/2):ycnt_indx2+ceil((no_glob_ygrid-1)/2)];
im2 = im2(indx_x2, indx_y2);
im2 = 10*log10(conv2(10.^(im2/10), gfilter,'same'));

im2_x = conv2(10.^(im2/10), Sx, 'same');
im2_y = conv2(10.^(im2/10), Sy, 'same');
im2_grad = 10*log10(sqrt(im2_x.^2+im2_y.^2));

%keyboard

%im1 = im1_grad;
%im2 = im2_grad;


% set parameters
nlevels = 5;
winsize = 3;
medfiltsize = 11;
nIterations = 5;

% use coarse to fine lucas kanade to obtain optical flow field
[u,v,warpI2] = coarse2fine_lk(im1,im2,nlevels,winsize,medfiltsize,nIterations);
                          
clear flow
flow(:,:,1) = u;
flow(:,:,2) = v;
                
offset = 30;
colorarray = 'wkkk';
%figure;imshow(im1);
%figure;imagesc(glob_x/1000, glob_y/1000, im1);
figure;imagesc(im1);
axis equal tight
caxis([40 120])
text(offset,offset,'Frame 1','FontSize',20,'FontWeight','bold','Color', 'k');
Frame(1) = getframe;
% for i = 1:-1:0
%     for j = 1:-1:0
%         % text(offset+j*2,offset+i*2,'Frame 1','FontSize',14,'FontName','Courier New','FontWeight','bold','Color',colorarray(i*2+j+1));
%         text(offset,offset,'Frame 1','FontSize',30,'FontWeight','bold','Color', 'k');
%     end
% end
% drawnow;
% Frame(1) = getframe;
% close;

figure;imagesc(im2);
axis equal tight
caxis([40 120])
text(offset,offset,'Frame 2','FontSize',20,'FontWeight','bold','Color', 'k');
Frame(2) = getframe;


%figure;imshow(warpI2);
%figure;imagesc(glob_x/1000, glob_y/1000, warpI2);
figure;imagesc(warpI2);
axis equal tight
caxis([40 120])
text(offset,offset,'Warped Frame 2','FontSize',20, 'FontWeight','bold','Color', 'k');
Frame(3) = getframe;
% for i = 1:-1:0
%     for j = 1:-1:0
%         % text(offset+j*2,offset+i*2,'Warped Frame 2','FontSize',14,'FontName','Courier New','FontWeight','bold','Color',colorarray(i*2+j+1));
%         text(offset,offset,'Warped Frame 2','FontSize',14,'FontName','Courier New','FontWeight','bold','Color', 'k');
%     end
% end
% drawnow;
% Frame(2) = getframe;
% close;

figure;
imshow(flowToColor(flow)); 
text(offset,offset,'Flow field','FontSize',20, 'FontWeight','bold','Color','k');
% for i = 1:-1:0
%     for j = 1:-1:0
%         % text(offset+j*2,offset+i*2,'Flow field','FontSize',20,'FontName','Courier New','FontWeight','bold','Color',colorarray(i*2+j+1));
%         text(offset,offset,'Flow field','FontSize',20,'FontName','Courier New','FontWeight','bold','Color','k');
%     end
% end
% drawnow;

flow_mag = sqrt(flow(:,:,1).^2+flow(:,:,2).^2);
figure; imagesc(flow_mag); colorbar



x = [-10:0.1:10];
y = x;
[X Y] = meshgrid(x,y);
uv_flow2color(:,:,1) = X;
uv_flow2color(:,:,2) = Y;
figure; imshow(flowToColor(uv_flow2color))
title('Flow to Color map')








%keyboard

fprintf('Press Ctrl+C to break the loop');
%h = figure;imagesc(glob_x/1000, glob_y/1000, Frame(1).cdata);
h = figure;imagesc(Frame(1).cdata);
axis equal tight
AxesH = get(h,'CurrentAxes');
for i = 1:100
    if mod(i,2) == 0
        %clf;
        cla(AxesH,'reset');
%        imshow(Frame(1).cdata,'Parent',AxesH);
%		imagesc(glob_x/1000, glob_y/1000, Frame(1).cdata,'Parent',AxesH);
		imagesc(Frame(1).cdata,'Parent',AxesH);
		axis equal tight
        drawnow;
    else
        %clf;
        cla(AxesH,'reset');
%        imshow(Frame(2).cdata,'Parent',AxesH);
%		imagesc(glob_x/1000, glob_y/1000, Frame(2).cdata,'Parent',AxesH);
		imagesc(Frame(3).cdata,'Parent',AxesH);
		axis equal tight
        drawnow;
    end
    pause(0.7);
end


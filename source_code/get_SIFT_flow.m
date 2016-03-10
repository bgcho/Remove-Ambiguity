function flow = get_SIFT_flow(im1, im2, SIFTflowpara)

im1=imresize(imfilter(im1,fspecial('gaussian',12,3),'same','replicate'),0.5,'bicubic');
im2=imresize(imfilter(im2,fspecial('gaussian',12,3),'same','replicate'),0.5,'bicubic');

im1=im2double(im1);
im2=im2double(im2);

%figure;imshow(im1);figure;imshow(im2);

cellsize=3;
gridspacing=1;

sift1 = mexDenseSIFT(im1,cellsize,gridspacing);
sift2 = mexDenseSIFT(im2,cellsize,gridspacing);

tic;[vx,vy,energylist]=SIFTflowc2f(sift1,sift2,SIFTflowpara);toc
vx = imresize(vx, 2, 'bicubic');
vy = imresize(vy, 2, 'bicubic');
flow(:,:,1)=vx;
flow(:,:,2)=vy;

%{
SIFTflowpara.alpha=0.01;
SIFTflowpara.d=1;
SIFTflowpara.gamma=0.001;
SIFTflowpara.nlevels=4;
SIFTflowpara.wsize=2;
SIFTflowpara.topwsize=10;
SIFTflowpara.nTopIterations = 60;
SIFTflowpara.nIterations= 40;
SIFTflowpara.wsize= 5;
%}

%       alpha:        (0.01) the weight of the truncated L1-norm regularization on the flow
%       d:            (1) the threshold of the truncation
%       gamma:        (0.001) the weight of the magnitude of the flow
%       nIterations:  (40) the number of iterations
%       nHierarchy:   (2)  the number of hierarchies in the efficient BP implementation
%       wsize:        (5)  the half size of the search window 




% warpI2=warpImage(im2,vx,vy);
% %figure;imshow(im1);figure;imshow(warpI2);
% figure;imagesc(im1);figure;imagesc(warpI2);
% title()
% xlabel()
% ylabel()
% caxis


% % display flow
% clear flow;


% figure;imshow(flowToColor(flow));
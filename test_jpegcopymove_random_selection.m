clear all;
close all;
addpath('jpegtbx_1.4');

%--------------%
%%% Settings %%%
imgPath = 'dataset\01_white_wall.jpg';  % the path of image file
maxNumOfIter = 50;	% max number of iteration
quality = 98;       % JPG quality
sizePatch = 16;     % size of testing patch
errnThresh = 1e-12; % threshold for stop iterative process
nComp = 3;          % number of component
numTest = 3;        % number of tests
ImgJPGName = 'imgtemp.jpg';  % temporary name of image for saving

%-------------------------%
%%% Read the image file %%%
Img = imread(imgPath);

%---------------------------------------------%
%%% Randmize position for SRC and DST patch %%%
imgBorder = 2 * sizePatch; % secure borders of image
srcPos = randi( [ imgBorder + 1, size( Img, 1 ) - imgBorder], 2, 1);
% position of src patch has to be aligned with the JPEG grid
srcPos = srcPos - mod( srcPos, 16 ) + 1;
% position of dst patch can be on a off-grid position
dstPos = randi( [ imgBorder + 1, size( Img, 1 ) - imgBorder], 2, numTest);

%----------------------------------------%
%%% Rewrite positions to struct format %%%
src = struct('x0',srcPos(1), 'y0',srcPos(2),'dx',sizePatch,'dy',sizePatch);
dst = struct('x0',dstPos(1,:), 'y0',dstPos(2,:),'dx',sizePatch,'dy',sizePatch);

%------------------------------%
%%% Changing the JPG quality %%%
imwrite(Img,ImgJPGName,'Quality',quality); % save image to new JPG quality

%-----------------------------%
%%% Call the main algorithm %%%
[result, numOfIter,L2,errn] = ...
    copymove_constraint(ImgJPGName, src, dst, maxNumOfIter, errnThresh,3);

%------------------------------%
%%% Visualization of results %%%
figure;
nRaws = floor(sqrt(numTest));
nCols = ceil(numTest/nRaws);

for i_test = 1:numTest     
    subplot(nRaws,nCols,i_test)
    plot(squeeze(errn(i_test,:,:))','LineWidth',3); 
    title(['Copy-move QCS, patch ' num2str(i_test)]);
    xlabel('Num. of Iterations');
    ylabel('Distance Between Sets');
    legend('Y'' component','Cb component','Cr component','Location','northeast');
end




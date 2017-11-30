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

%---------------------------------------%
%%% Randmize position for SRC and DST %%%
imgBorder = 2 * sizePatch; % secure borders of image
srcPos = randi( [ imgBorder + 1, size( Img, 1 ) - imgBorder], 2, 1);
% position of src patch has to be aligned with the JPEG grid
srcPos = srcPos - mod( srcPos, 16 ) + 1;
% rewrite SRC positions to struct format
src = struct('x0',srcPos(1), 'y0',srcPos(2),'dx',sizePatch,'dy',sizePatch);

%%% Randmize position for copy-move DST patch %%%
% position of dst patch can be on a off-grid position
possibleDstRow(1,:) = [ imgBorder + 1 : srcPos(1) - sizePatch - 1, ...
    sizePatch +  srcPos(1) + 1 :  size( Img, 1 ) - imgBorder];
possibleDstCol(2,:) = [ imgBorder + 1 : srcPos(2) - sizePatch - 1, ...
    sizePatch +  srcPos(2) + 1 :  size( Img, 1 ) - imgBorder];
% control that we do not rewrite the DST parch to SRC patch
dstPos = zeros(2,numTest);
for i_test = 1 : numTest
    dstPos(1,i_test) = possibleDstRow(randi(size(possibleDstRow,2)));
    dstPos(2,i_test) = possibleDstRow(randi(size(possibleDstRow,2)));
end

%-------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% COPY-MOVE TAMPERING %%%%%%
for i_test = 1 : numTest    
    tempImg = Img;
    dst = struct('x0',dstPos(1,i_test), 'y0',dstPos(2,i_test),'dx',sizePatch,'dy',sizePatch);
    cropImg = Img(src.y0:src.y0+src.dy-1,src.x0:src.x0+src.dx-1,:);
    tempImg(dst.y0:dst.y0+dst.dy-1,dst.x0:dst.x0+dst.dx-1,:) = cropImg;
    imwrite(tempImg,ImgJPGName,'Quality',quality); % save image to new JPG quality
    %%% Call the main algorithm %%%
    [result(i_test), numOfIter(i_test,:), L2(i_test), errn(i_test,:,:)] = ...
        copymove_constraint(ImgJPGName, src, dst, maxNumOfIter, errnThresh,3);
end
%%% END OF COPY-MOVE TAMPERING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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




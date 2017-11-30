function [result, numIter, L2, errn] = copymove_constraint(fname,src,dst,nIters,errnThresh,nComp)
%COPYMOVE_CONSTRAINT Tests possibility of copy-move forgery in JPEG files.
%   [result, numIter, L2, errn] = copymove_constraint(fname,src,dst,nIters,errnThresh,nComp)
%   INPUTS:
%   fname       ...name of test image
%               fname = 'image.jpg'
%   src         ...position of source patch
%               src = struct('x0', ROWPOS, 'y0', COLPOS, 'dx', SIZEPATCH, 'dy', SIZEPATCH);
%   dst         ...position of destination patch
%               dst = struct('x0', ROWPOS, 'y0', COLPOS, 'dx', SIZEPATCH, 'dy', SIZEPATCH);
%   nIter       ...the max number of iteration
%               nIter = 50 (default)
%   errnThresh  ...threhold for minimum distance between sets
%               errnThresh = 50 (default)
%   nComp       ...the number of components we want to compute
%               nComp = 1 ... only the Y' component
%               nComp = 3 ... for all three components - Y'CbCr (default) 
%   OUTPUTS:
%   result      the end distance between sets for Y' component
%   numIter     how muny iteration we need for reached the errnThresh
%   L2          L2 distance between src patch and dst patches
%   errn        the distance for required components in each iteration
%   
%   Michal Sorel & Adam Novozamsky 2016 

    % Set up default parameters
    if ~exist('nComp','var'), nComp = 1; end
    if ~exist('nIters','var'), nIters = 50; end
    if ~exist('errnThresh','var'), errnThresh = 1e-12; end
    
    if nIters < 1
        error('The number of iteration must be greater than 0.');
    end
    
    if errnThresh < 0 && errnThresh > 1
        error('The threshold ??? must be in range (0,1).');
    end

    kernel=1;
    aux = imread(fname);   
    y = jpeg_read(fname);
    
    if nComp > y.image_components
        error(['The required number of components exceeds the number of components in the image. '...
            'In this case the max number of component = ' num2str(y.image_components)]);
    end
       
    K = cell(y.image_components,1);
    Q = cell(y.image_components,1);
    k = cell(y.image_components,1);
    lower_bound_src = cell(y.image_components,1);
    upper_bound_src = cell(y.image_components,1);
    lower_bound_dest = cell(y.image_components,1);
    upper_bound_dest = cell(y.image_components,1);
    K_src = cell(y.image_components,1); 
    K_dest = cell(y.image_components,1);
    size_downsampled_padded = cell(y.image_components,1);
    sampling_factor = zeros(y.image_components,2);
    factor = zeros(y.image_components,2);
    numIter = zeros(length(dst.x0),y.image_components);
    result = zeros(1,length(dst.x0));
    errn = zeros(length(dst.x0),y.image_components, nIters);
    L2 = zeros(1,length(dst.x0));
    for i_comp = 1 : y.image_components
            sampling_factor(i_comp,:)=[y.comp_info(i_comp).v_samp_factor y.comp_info(i_comp).h_samp_factor];
            factor(i_comp,:)=sampling_factor(i_comp,:)./max(sampling_factor);
            if mod(src.x0,8/factor(i_comp,2))~=1 || mod(src.y0,8/factor(i_comp,1))~=1
                error('Source rectangle must be aligned with 8x8 grid (16x16 for color).'); 
            end             
            if kernel
                k{i_comp}=factor(i_comp,1)*factor(i_comp,2);
            end    
            size_downsampled_padded{i_comp}=size(y.coef_arrays{i_comp}); % get size of downsapled channels

            % replicated inverse quantization table over whole image
            Q{i_comp}=repmat(1./y.quant_tables{y.comp_info(i_comp).quant_tbl_no},...
                size_downsampled_padded{i_comp}/size(y.quant_tables{y.comp_info(i_comp).quant_tbl_no},1)); 
            K{i_comp} = 1./(k{i_comp}*Q{i_comp}.^2);
    end
    x_est = cell(y.image_components,1);
    dest = src;
    for i_dest = 1:length(dst.x0)
        dest.x0 = dst.x0(i_dest);
        dest.y0 = dst.y0(i_dest);
        L2(i_dest) = norm(vec(double(getcrop(aux,src)-getcrop(aux,dest))));
        %disp(['Initial ||src-dest||= ' num2str(norm(double(getcrop(aux,src)-getcrop(aux,dest)),'fro'))])
%         disp(['RGB initial ||src-dest||= ' num2str(norm(vec(double(getcrop(aux,src)-getcrop(aux,dest))),'fro'))])       
%         disp('--------------------------------------------------');
        d = cell(size(x_est));
        a = cell(size(x_est));
       % niters = 12;        
        err = zeros(y.image_components,1);
%         errn = zeros(size(err));        
        
        z_ycbcr = zeros(y.image_height,y.image_width,y.image_components);             
        for i_comp = 1 : y.image_components 
            maxDiff = zeros(1,3);
            destL = destEnlarge(dest,factor(i_comp,:));    
            srcdown = rectdown(src,factor(i_comp,:));
            destLdown = rectdown(destL,factor(i_comp,:));
            crop = struct('x0',dest.x0-destL.x0+1,'y0',dest.y0-destL.y0+1,'dx',src.dx,'dy',src.dy);                    
            lower_bound_src{i_comp} = getcrop(y.coef_arrays{i_comp},srcdown)-1/2;
            upper_bound_src{i_comp} = getcrop(y.coef_arrays{i_comp},srcdown)+1/2-100*eps; % to treat open interval
            lower_bound_dest{i_comp} = getcrop(y.coef_arrays{i_comp},destLdown)-1/2;
            upper_bound_dest{i_comp} = getcrop(y.coef_arrays{i_comp},destLdown)+1/2-100*eps; % to treat open interval
            K_src{i_comp} = getcrop(K{i_comp},srcdown);
            K_dest{i_comp} = getcrop(K{i_comp},destLdown);             
            % standard decompresion
            aux=ibdct(dequantize(y.coef_arrays{i_comp},y.quant_tables{y.comp_info(i_comp).quant_tbl_no}));
            aux = resample_up_down(aux,size(aux)./factor(i_comp,:),kernel)/prod(factor(i_comp,:));
            z_ycbcr(:,:,i_comp) = aux(1:y.image_height,1:y.image_width);             
            
            err(i_comp) = norm(double(getcrop(z_ycbcr(:,:,i_comp),src)-getcrop(z_ycbcr(:,:,i_comp),dest)),'fro');
%             disp(['Channel ' num2str(n) ' initial ||src-dest||= ' num2str(err(n))]);            
            if err(i_comp) >= errnThresh               
                % initialize ADMM auxiliary variables
                a{i_comp} = double(getcrop(z_ycbcr(:,:,i_comp),destL));
                d{i_comp} = zeros(size(a{i_comp})); 
                for i_iter = 1:nIters 
                    numIter(i_dest,i_comp) = numIter(i_dest,i_comp) + 1;
                    x_est{i_comp} = proj2QCSmasked(a{i_comp}+d{i_comp},[],K_dest{i_comp},getcrop(Q{i_comp},destLdown),[destLdown.dy destLdown.dx],...
                                                            lower_bound_dest{i_comp}, upper_bound_dest{i_comp},kernel);
                    a{i_comp} = proj2QCSmasked(x_est{i_comp}-d{i_comp},crop,K_src{i_comp},getcrop(Q{i_comp},srcdown),[srcdown.dy srcdown.dx],...
                                                            lower_bound_src{i_comp}, upper_bound_src{i_comp},kernel);
                    d{i_comp} = d{i_comp} - (x_est{i_comp} - a{i_comp});                                         
                    errn(i_dest,i_comp,i_iter) = norm(x_est{i_comp}-a{i_comp},'fro');                   
                    if errn(i_dest,i_comp,i_iter) < errnThresh, break; end                                                                                        
                end
            else
                errn(i_dest,i_comp, 1)=0;
            end
        end
        x_res = zeros(size(x_est{1},1),size(x_est{1},2),y.image_components);
        a_res = zeros(size(x_res));       
        for i_comp = 1:y.image_components
            if err(i_comp) >=errnThresh
                maxDiff(i_comp) = max(abs(x_est{i_comp}(:)-a{i_comp}(:)));                
            else
                maxDiff(i_comp) = 0;                                
            end
            result(i_dest) = result(i_dest) + errn(i_dest,i_comp,max(1,numIter(i_dest,i_comp)));
        end 
    end
    errn = squeeze(errn);
end 

function r = rectdown(sr,factor) % rectangle in down-sample coordinates
    r = sr;
    r.x0 = (sr.x0-1)*factor(2)+1;
    r.y0 = (sr.y0-1)*factor(1)+1;
    r.dx = sr.dx*factor(2);
    r.dy = sr.dy*factor(1);
end

function destL = destEnlarge(dest,fact)
%destEnlarge Enlarges the destination are structure to nearest multiple of 8x8 pixels
%used in JPEGCopyMoveTest

    destL = dest; % large dest
    f = 8./fact;
    destL.x0 = dest.x0-mod(dest.x0-1,f(2));
    destL.y0 = dest.y0-mod(dest.y0-1,f(1));
    if mod(dest.x0-1,f(2)) == 0
        destL.dx = dest.dx;
    else
        destL.dx = dest.dx + f(2);
    end
    if mod(dest.y0-1,f(1)) == 0
        destL.dy = dest.dy;
    else
        destL.dy = dest.dy + f(1);
    end
end

function r = getcrop(x,src)
    r = x(src.y0:src.y0+src.dy-1,src.x0:src.x0+src.dx-1,:);
end

function V = vec(A)
% M ... matrix MxN
% V ... vector (M*N)x1
    V = A(:);
end

function y=resample_up_down(x,size_y,kernel)
% simple resampling up/down for JPEG purposes
    y=zeros(size_y);
    switch kernel
        case 0  % pure resample
            if size(x)<size_y; % up

                y(1:2:end,1:2:end)=x;
            elseif size(x)>size_y % down
                y=x(1:2:end,1:2:end);
            else
                y=x;
            end

        case 1 % with filter
            h=[1,1,0;1,1,0;0 0 0]/4;
    %         h=imrotate(h,180);        
            if size(x)>size_y; % down

    %             y=[[y(1,1) y(1,:)];[y(:,1) y]];
                y=[[y y(:,end)];[ y(end,:) y(end,end)]];
                y=convn(x,h,'same');
                y=y(1:2:end,1:2:end);
    %             y=y(1:2:end,1:2:end)/4;            
            elseif size(x)<size_y % up
                h=imrotate(h,180);
                y(1:2:end,1:2:end)=x;
    %             
                y=convn(y,h,'same');
            else
                y=x;
            end
    end
end

function x_est = proj2QCSmasked(y,crop,K,Q,size_downsampled,lower_bound,upper_bound,kernel)            
%
%Call: x_est{m} = proj2QCSmasked(y,[],K{m},Q{m},size_downsampled{m},size_downsampled_padded{m}, lower_bound{m}, upper_bound{m})
%    
    if isempty(crop)
        rs = resample_up_down(y,size_downsampled,kernel); % downsample
    else
        rs = resample_up_down(getcrop(y,crop),size_downsampled,kernel);                                                
    end
    rs = Q.*bdct(rs);
    rs = rs - box_projection(rs, lower_bound, upper_bound);
    rs = K.*rs;
    rs = ibdct(rs.*Q);
    %rs = rs(1:size_downsampled(1),1:size_downsampled(2)); % crop padded values   
    if isempty(crop)
        rs = resample_up_down(rs,size(y),kernel); % upsample
        x_est = y - rs;
    else
        rs = resample_up_down(rs,[crop.dy crop.dx],kernel); % upsample
        aux = zeros(size(y));
        aux(crop.y0:crop.y0+crop.dy-1,crop.x0:crop.x0+crop.dx-1) = rs;
        x_est = y - aux; 
    end
end

function x=box_projection(x,lower_bound,upper_bound)
% Projection of 2D data x onto its constraints defined by lower and upper bounds
% For projections on quantization constraint sets
x = min(cat(3, max(cat(3,lower_bound,x),[],3), upper_bound),[],3);
end
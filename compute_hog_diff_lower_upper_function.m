function [upperhogdiff,lowerhogdiff]=compute_hog_diff_lower_upper_function(upperatlas,loweratlas, ...
        cellsize,endslice,numexperimental,indices,experimental_file,startslice,tatlas,tmask,experimental_thickness,...
thick_thresh,maskfile,expresolution)

upperhogdiff=zeros(numexperimental,size(upperatlas,3));
lowerhogdiff=zeros(numexperimental,size(loweratlas,3));
for i=1:numexperimental
    i
    experimental_index=indices(1,i);
    exp_slice=single(experimental_file(:,:,experimental_index));
    small_exp=imresize(exp_slice,expresolution/25);
    mask_slice=logical(maskfile(:,:,experimental_index));
    small_mask=imresize(mask_slice,expresolution/25,'nearest');
    startposition=startslice+(experimental_index-indices(1,1))*experimental_thickness/25*(1-thick_thresh);
    endposition=endslice-(indices(1,end)-experimental_index)*experimental_thickness/25*(1-thick_thresh);
    parfor k=round(startposition):min(endslice,round(endposition))
        %warp before cut
        atlasslice=tatlas(:,:,k);
        %warp to downsized experimental slice
        atlasmask = imfill(atlasslice>0,'holes');
        atlasmask = bwareafilt(atlasmask,1);
        [~,tform]=warp_with_shape(atlasmask,bwareafilt(small_mask,1));
        %warp the mask of the atlas file
        warpedatlas = imtransform(atlasslice,tform,'bilinear','XData',[1 size(small_mask,2)],'YData',...
            [1 size(small_mask,1)]);
        warpedmask=imtransform(tmask(:,:,k),tform,'nearest','XData',[1 size(warpedatlas,2)],'YData',[1 size(warpedatlas,1)]);
        %cut with the original cut
        warpedatlasupper=warpedmask.*warpedatlas;
        warpedatlaslower=(~warpedmask).*warpedatlas;%(~tmask(:,:,k)).*warpedexp;
        upperatlashog=vl_hog(single(warpedatlasupper),cellsize);
        loweratlashog=vl_hog(single(warpedatlaslower),cellsize);
        [m,n,l]=size(upperatlashog);
        upperhogexp=vl_hog(single(small_exp.*warpedmask),cellsize);
        lowerhogexp=vl_hog(single(small_exp.*(~warpedmask)),cellsize);
		diffupper=upperhogexp-upperatlashog;
		upperhogdiff(i,k)=norm(reshape(diffupper,1,m*n*l));
		difflower=lowerhogexp-loweratlashog;
		lowerhogdiff(i,k)=norm(reshape(difflower,1,m*n*l));
     end
end
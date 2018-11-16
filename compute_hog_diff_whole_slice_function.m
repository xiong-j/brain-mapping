function hogdiff=compute_hog_diff_whole_slice_function(tatlas, cellsize,...
    endslice,numexperimental,origindices,experimental_file, startslice,...
mask_file, matching_min, matching_max,expResolution)

hogdiff=zeros(numexperimental,size(tatlas,3));
for i=1:numexperimental
    i
    experimental_index=origindices(1,i);
    exp_slice=single(experimental_file(:,:,experimental_index))/255;
    small_exp=imresize(exp_slice,expResolution/25);% to the resolution of atlas
    mask_slice=double(mask_file(:,:,experimental_index))/255;
    small_mask=imresize(mask_slice,expResolution/25);
    for k = max(matching_min(1,origindices(i)), startslice): min(matching_max(1,origindices(i)), endslice)
        atlasslice=tatlas(:,:,k);
        %warp to downsized experimental slice
        atlasmask = imfill(atlasslice>0,'holes');
        [warpedatlas,tform] = warp_with_shape(atlasmask,small_mask);
        warpedatlas=imtransform(atlasslice,tform,'XData',[1 size(small_mask,2)],'YData',[1 size(small_mask,1)]);
        atlashog=vl_hog(single(warpedatlas),cellsize);
        [m,n,l]=size(atlashog);
        hogexp=vl_hog(single(small_exp),cellsize);
        mask = ones(m,n,l);
        mask(round(m/4):m,:,l)=2;
		this_diff=(hogexp-atlashog).*mask;
		hogdiff(i,k)=norm(reshape(this_diff,1,m*n*l));
     end
end
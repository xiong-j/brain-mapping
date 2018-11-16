function [output_file,mat_file,x,y,correspondences]=find_slice_matches_function(experimental_thickness, ...
    experimental_file, x,y,maskfile, expResolution,brainname,atlas_folder,...
    good_quality_indices,matching1,matching2,startslice,endslice,thresh)
output_file = strcat(brainname,'_correspondences','.tif');
mat_file = strcat(brainname,'_correspondences','.mat');
neighbor_size = 10;
cellsize=15;
experimental_volume_size=size(experimental_file,3);
if y>=0
    filename=strcat(atlas_folder,'/tatlas_',int2str(x),'_',int2str(y),'.mat');
else
    filename=strcat(atlas_folder,'/tatlas_',int2str(x),'_',int2str(-y),'_negative.mat');
end
load([filename])
matching_min = zeros(1,experimental_volume_size);% all the experimental slices
matching_max = matching_min;
num_experimental = length(good_quality_indices);
matching_min(1,good_quality_indices)=min(matching1,matching2)-neighbor_size;
matching_max(1,good_quality_indices)=max(matching1,matching2)+neighbor_size;
% check if starting slice is 1
if min(good_quality_indices)~=1
    matching_max(1,1) = round(matching_max(1,min(good_quality_indices))-experimental_thickness/25*(min(good_quality_indices)-1)*(1-thresh)+neighbor_size);
    matching_min(1,1) = round(matching_min(1,min(good_quality_indices))-experimental_thickness/25*(min(good_quality_indices)-1)*(1+thresh)-neighbor_size);
    good_quality_indices = [1 good_quality_indices];
    num_experimental = num_experimental+1;
end
for i = 1:length(experimental_volume_size)
    if matching_min(1, i) == 0
        matching_min(1,i) = matching_min(1,i-1);
    end
end
 %find the last nonzero
if max(good_quality_indices)~=size(experimental_file,3)
    matching_max(1,end)=round(matching_max(1,max(good_quality_indices))+experimental_thickness/25*(size(experimental_file,3)...
        -max(good_quality_indices))*(1+thresh)+neighbor_size);
    matching_min(1,end)=round(matching_min(1,max(good_quality_indices))+experimental_thickness/25*(size(experimental_file,3)...
        -max(good_quality_indices))*(1-thresh)-neighbor_size);
    good_quality_indices=[good_quality_indices,size(experimental_file,3)];
    num_experimental=num_experimental+1;
end

hogdiff=compute_hog_diff_whole_slice_function(tatlas, cellsize,endslice, ...
    num_experimental , good_quality_indices,experimental_file, startslice,...
    maskfile, matching_min, matching_max,expResolution);
[sumhog,answer]=dp_hog_function(endslice,hogdiff(:,1:endslice),num_experimental,good_quality_indices,experimental_thickness);
correspondences=round(interp1(good_quality_indices,answer,1:size(experimental_file,3)));
for i=1:length(correspondences)
    resultimg=tatlas(:,:,correspondences(i));
    expimg=experimental_file(:,:,i)/255;
    expimg=imresize(expimg,size(resultimg));
    imwrite(imadjust(expimg),output_file,'WriteMode','append');
    imwrite(imadjust(resultimg),output_file,'WriteMode','append');
end
indices = good_quality_indices;
save(mat_file,'answer','indices','correspondences','y','x');
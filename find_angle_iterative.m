function [y,x,matching1,matching2,startslice,endslice]=find_angle_iterative(experimental_thickness,maskfile,brainfile,...
    begin_angle,indices,atlas_folder,expresolution)
% Use the area of the first and last slice of the input file to estimate
% the searching range

num_experimental=size(maskfile,3);
thick_thresh=0.2;
experimental_file=brainfile;
numexperimental=length(indices);
cellsize=15;
y=begin_angle;
while 1
    tic
    x=0;
    y
    if y>=0
        filename=strcat(atlas_folder,'/','tatlas_0_',int2str(y),'.mat');
    else
        filename=strcat(atlas_folder,'/','tatlas_0_',int2str(-y),'_negative.mat');
    end
    %check if the rotated file exist
    if exist(filename, 'file') == 2
        load([filename])
    else
        load(strcat(atlas_folder,'/tatlas_0_0.mat'))
        [tatlas,origin]=rotate_w_label(tatlas,x,y);
        save(strcat(filename),'tatlas','origin','-v7.3');
    end
    tmask=zeros(size(tatlas));
    for i=1:size(tatlas,3);
        [r,c]=find(tatlas(:,:,i));
        upperline=min(r);
        lowerline=max(r);
        midline=round((upperline+lowerline)/2);
        tmask(1:midline-1,:,i)=1;
    end
    upperatlas=tatlas.*tmask;
    loweratlas=tatlas.*(1-tmask);
    first_slice=single(experimental_file(:,:,1));
    %resize to the atlas resolution
    resizedFirstSlice = imresize(first_slice,expresolution/25);
    firstHeight = length(find(sum(resizedFirstSlice,2)));
    %% use centroid instead
    flag = 1;
    startslice = 0;
    while flag == 1
        startslice = startslice+1;
        atlasslice = tatlas(:,:,startslice);
        atlasslice=bwareafilt(atlasslice>0,1).*atlasslice;
        stats = regionprops(atlasslice>0,'BoundingBox');
        if ~isempty(stats)
            bb = stats.BoundingBox;
            heightstart = bb(2);
            heightend = bb(2)+bb(4);
            thisHeight = heightend-heightstart+1;
            if thisHeight>(1-thick_thresh)*firstHeight && thisHeight<(1+thick_thresh)*firstHeight
                heightmid = (heightstart+heightend)/2;
                %if heightmid>=0.4*size(tatlas,1) && heightmid<=0.6*(size(tatlas,1)) && heightend-heightstart>50
                    flag=0;
                %end
            end
        end
    end
    flag = 1;
    endslice = size(tatlas,3)+1;
    last_slice=single(experimental_file(:,:,num_experimental));
    resizedLastSlice = imresize(last_slice,expresolution/25);
    lastHeight = length(find(sum(resizedLastSlice,2)));
    while flag == 1
        endslice = endslice-1;
        atlasslice = tatlas(:,:,endslice);
        atlasslice=bwareafilt(atlasslice>0,1).*atlasslice;
        stats = regionprops(atlasslice>0,'BoundingBox');
        if ~ isempty(stats)
            bb=stats.BoundingBox;
            heightstart = bb(2);
            heightend = bb(2)+bb(4);
            heightmid = (heightstart+heightend)/2;
            thisHeight = heightend-heightstart;
            if thisHeight>(1-thick_thresh)*lastHeight && thisHeight<(1+thick_thresh)*lastHeight
                %if heightmid>=0.4*size(tatlas,1) && heightmid<=0.6*(size(tatlas,1)) && heightend-heightstart>50
                    flag = 0;
                %end
            end
        end
    end
    [upperhogdiff,lowerhogdiff]=compute_hog_diff_lower_upper_function(upperatlas,loweratlas, ...
        cellsize,endslice,numexperimental,indices,experimental_file, ...
        startslice,tatlas,tmask,experimental_thickness,...
        thick_thresh,maskfile,expresolution);
    [~,upperanswer]=dp_hog_function(endslice,upperhogdiff(:,1:endslice),numexperimental,indices,experimental_thickness);
    [~,loweranswer]=dp_hog_function(endslice,lowerhogdiff(:,1:endslice),numexperimental,indices,experimental_thickness);
    thisdiff = mean(upperanswer-loweranswer);
    if y==begin_angle
        lastsign = sign(thisdiff);
    else
        if lastsign~=sign(thisdiff)
            break
        end
    end
    lasty=y;
    if thisdiff>0
        y=y-1;
    else
        y=y+1;
    end
    lastdiff=thisdiff;
    toc
end
if abs(lastdiff)<abs(thisdiff)
    y=lasty;
end
x=0;
while 1
    tic
    x
    if y>=0
        filename=strcat(atlas_folder,'/','tatlas_',int2str(x),'_',int2str(y),'.mat');
    else
        filename=strcat(atlas_folder,'/','tatlas_',int2str(x),'_',int2str(-y),'_negative.mat');
    end
    %check if the rotated file exist
    if exist(filename, 'file') == 2
        load([filename])
    else
        load(strcat(atlas_folder,'/tatlas_0_0.mat'))
        [tatlas,origin]=rotate_w_label(tatlas,x,y);
        save(strcat(filename),'tatlas','origin','-v7.3');
    end
    tmask=zeros(size(tatlas));
    for i=1:size(tatlas,3);
        [r,c]=find(tatlas(:,:,i));
        rightline = max(c);
        leftline = min(c);
        midline=round((leftline+rightline)/2);
        tmask(:,1:midline-1,i)=1;%mask the left
    end
    leftatlas=tatlas.*tmask;
    rightatlas=tatlas.*(1-tmask);
    first_slice=single(experimental_file(:,:,1));
    %resize to the atlas resolution
    resizedFirstSlice = imresize(first_slice,expresolution/25);
    firstHeight = length(find(sum(resizedFirstSlice,2)));
    %% use centroid instead
    flag = 1;
    startslice = 0;
    while flag == 1
        startslice = startslice+1;
        atlasslice = tatlas(:,:,startslice);
        atlasslice = bwareafilt(atlasslice>0,1).*atlasslice;
        stats = regionprops(atlasslice>0,'BoundingBox');
        if ~isempty(stats)
            bb = stats.BoundingBox;
            heightstart = bb(2);
            heightend = bb(2)+bb(4);
            thisHeight = heightend-heightstart+1;
            if thisHeight>(1-thick_thresh)*firstHeight && thisHeight<(1+thick_thresh)*firstHeight
                heightmid = (heightstart+heightend)/2;
                %if heightmid>=0.4*size(tatlas,1) && heightmid<=0.6*(size(tatlas,1)) && heightend-heightstart>50
                    flag=0;
                %end
            end
        end
    end
    flag = 1;
    endslice = size(tatlas,3)+1;
    last_slice=single(experimental_file(:,:,num_experimental));
    resizedLastSlice = imresize(last_slice,expresolution/25);
    lastHeight = length(find(sum(resizedLastSlice,2)));
    while flag == 1
        endslice = endslice-1;
        atlasslice = tatlas(:,:,endslice);
        atlasslice = bwareafilt(atlasslice>0,1).*atlasslice;
        stats = regionprops(atlasslice>0,'BoundingBox');
        if ~ isempty(stats)
            bb=stats.BoundingBox;
            heightstart = bb(2);
            heightend = bb(2)+bb(4);
            heightmid = (heightstart+heightend)/2;
            thisHeight = heightend-heightstart;
            if thisHeight>(1-thick_thresh)*lastHeight && thisHeight<(1+thick_thresh)*lastHeight
                %if heightmid>=0.4*size(tatlas,1) && heightmid<=0.6*(size(tatlas,1)) && heightend-heightstart>50
                    flag = 0;
                %end
            end
        end
    end
    [upperhogdiff,lowerhogdiff]=compute_hog_diff_lower_upper_function(upperatlas,loweratlas, ...
        cellsize,endslice,numexperimental,indices,experimental_file, ...
        startslice,tatlas,tmask,experimental_thickness,...
        thick_thresh,maskfile,expresolution);
    [~,matching1]=dp_hog_function(endslice,upperhogdiff(:,1:endslice),numexperimental,indices,experimental_thickness);
    [~,matching2]=dp_hog_function(endslice,lowerhogdiff(:,1:endslice),numexperimental,indices,experimental_thickness);
    thisdiff = mean(matching1-matching2);
    if x==0
        lastsign = sign(thisdiff);
    else
        if lastsign~=sign(thisdiff)
            break
        end
    end
    lastx=x;
    if thisdiff>0
        x=x+1;
    else
        x=x-1;
    end
    lastdiff=thisdiff;
    laststartslice = startslice;
    lastendslice=endslice;
    lastmatching1 = matching1;
    lastmatching2 = matching2;
    toc
end
if abs(lastdiff)<abs(thisdiff)
    x=lastx;
    startslice=laststartslice;
    endslice = lastendslice;
    matching1 = lastmatching1;
    matching2 = lastmatching2;
end
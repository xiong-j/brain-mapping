function check_angle_result(outputfile, threshold, brainfile, angle_yrange)
load([outputfile]);
tmp=diffs(2:end,:);
avg_diffs=abs(mean(tmp,1));
best_index=find(avg_diffs<threshold);%here is just the index
minangle=min(angle_yrange);
for i=best_index
    y=angle_yrange(i);
    % answertop starts from the min angle
    answertop=answerstop(1,y-minangle+1);
    answerbottom=answersbottom(1,y-minangle+1);
    if y>=0
        filename=strcat('tatlas_0_',int2str(y),'.mat');
    else
        filename=strcat('tatlas_0_',int2str(-y),'_negative.mat');
    end
    load([filename])
    answertop=answertop{1};
    answerbottom=answerbottom{1};
    l=size(answertop,1);
    [pathstr,name,ext]=fileparts(outputfile);
    [m,n,slicenum]=size(tatlas);
    for k=1:l
            img=double(brainfile(:,:,indices(k)));
            img=imresize(img,[m,n]);
            result_name=strcat(name,'half_matching_result',int2str(y),'.tif');
            imwrite(imadjust(img),result_name,'WriteMode','append');
            imwrite(imadjust(tatlas(:,:,answertop(k))),result_name,'WriteMode','append');
            imwrite(imadjust(tatlas(:,:,answerbottom(k))),result_name,'WriteMode','append');
    end
end
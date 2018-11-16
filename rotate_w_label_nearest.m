function [rotated,origin]=rotate_w_label_nearest(volume,x,y)
% rotate the atlas and keep tracks of where each pixel comes from

thx=x*pi/180;
thy=y*pi/180;
R=makehgtform('xrotate',thx,'yrotate',thy);
[rotated,origin]=find_origin_nearest(volume,R);

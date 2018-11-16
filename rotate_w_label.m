%check find origin
function [rotated,origin]=rotate_w_label(volume,x,y)
thx=x*pi/180;
thy=y*pi/180;
R=makehgtform('xrotate',thx,'yrotate',thy);
[rotated,origin]=find_origin(volume,R);

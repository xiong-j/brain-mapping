function [new_img, origin] = find_origin_nearest(old_img, old_M, new_elem_size, verbose, bg, method)
% rotate a volume and with the origin of each pixel

   if ~exist('old_img','var') | ~exist('old_M','var')
      error('Usage: [new_img new_M] = affine(old_img, old_M, [new_elem_size], [verbose], [bg], [method]);');
   end

   if ndims(old_img) == 3
      if ~isequal(size(old_M),[4 4])
         error('old_M should be a 4x4 affine matrix for 3D volume.');
      end
   elseif ndims(old_img) == 2
      if ~isequal(size(old_M),[3 3])
         error('old_M should be a 3x3 affine matrix for 2D image.');
      end
   else
      error('old_img should be either 2D image or 3D volume.');
   end

   if ~exist('new_elem_size','var') | isempty(new_elem_size)
      new_elem_size = [1 1 1];
   elseif length(new_elem_size) < 2
      new_elem_size = new_elem_size(1)*ones(1,3);
   elseif length(new_elem_size) < 3
      new_elem_size = [new_elem_size(:); 1]';
   end

   if ~exist('method','var') | isempty(method)
      method = 1;
   elseif ~exist('bresenham_line3d.m','file') & method == 3
      error([char(10) char(10) 'Please download 3D Bresenham''s line generation program from:' char(10) char(10) 'http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=21057' char(10) char(10) 'to test Fischer''s Bresenham interpolation method.' char(10) char(10)]);
   end

   %  Make compatible to MATLAB earlier than version 7 (R14), which
   %  can only perform arithmetic on double data type
   %
   old_img = double(old_img);
   old_dim = size(old_img);

   if ~exist('bg','var') | isempty(bg)
      bg = mean([old_img(1) old_img(end)]);
   end

   if ~exist('verbose','var') | isempty(verbose)
      verbose = 1;
   end

   if ndims(old_img) == 2
      old_dim(3) = 1;
      old_M = old_M(:, [1 2 3 3]);
      old_M = old_M([1 2 3 3], :);
      old_M(3,:) = [0 0 1 0];
      old_M(:,3) = [0 0 1 0]';
   end

   %  Vertices of img in voxel
   % make original vertices
   XYZvox = [	1		1		1
		1		1		old_dim(3)
		1		old_dim(2)	1
		1		old_dim(2)	old_dim(3)
		old_dim(1)	1		1
		old_dim(1)	1		old_dim(3)
		old_dim(1)	old_dim(2)	1
		old_dim(1)	old_dim(2)	old_dim(3)   ]';
%rotation
   old_R = old_M(1:3,1:3);
   old_T = old_M(1:3,4);
%translation
   %  Vertices of img in millimeter
   % vertices after transformation
   %(3x3)x(3x8)=(3x8)+translation, warp vertices
   XYZmm = old_R*(XYZvox-1) + repmat(old_T, [1, 8]);

   %  Make scale of new_M according to new_elem_size
   % new m now diagonal matrix 4x4
   new_M = diag([new_elem_size 1]);

   %  Make translation so minimum vertex is moved to [1,1,1]
   %XYZmm vertices after transformation
   new_M(1:3,4) = round( min(XYZmm,[],2) );
%take minimum along rows.

   %  New dimensions will be the maximum vertices in XYZ direction (dim_vox)
   %  i.e. compute   dim_vox   via   dim_mm = R*(dim_vox-1)+T
   %  where, dim_mm = round(max(XYZmm,[],2));
   %
   new_dim = ceil(new_M(1:3,1:3) \ ( round(max(XYZmm,[],2))-new_M(1:3,4) )+1)';
%new max vertex
%new M is just a diagonal matrix
%R is now just no rotation 
%old dim is size of old img
%new dim is size.all rough

   %  Initialize new_img with new_dim
   %
   new_img = zeros(new_dim(1:3));
   origin=zeros([new_dim(1:3) 3]);
   

   %  Mask out any changes from Z axis of transformed volume, since we
   %  will traverse it voxel by voxel below. We will only apply unit
   %  increment of mask_Z(3,4) to simulate the cursor movement
   %
   %  i.e. we will use   mask_Z * new_XYZvox   to replace   new_XYZvox
   %
   mask_Z = diag(ones(1,4));%diagonal identity matrix 4x4
   mask_Z(3,3) = 0;%mask out z change

   %  It will be easier to do the interpolation if we invert the process
   %  by not traversing the original volume. Instead, we traverse the
   %  transformed volume, and backproject each voxel in the transformed 
   %  volume back into the original volume. If the backprojected voxel
   %  in original volume is within its boundary, the intensity of that
   %  voxel can be used by the cursor location in the transformed volume.
   %
   %  First, we traverse along Z axis of transformed volume voxel by voxel
   %
   for z = 1:new_dim(3)

      if verbose & ~mod(z,10)
         fprintf('%.2f percent is done.\n', 100*z/new_dim(3));
      end

      %  We need to find out the mapping from voxel in the transformed
      %  volume (new_XYZvox) to voxel in the original volume (old_XYZvox)
      %XYZmm is eight vertices after transformation
      %  The following equation works, because they all equal to XYZmm:
      %  new_R*(new_XYZvox-1) + new_T  ==  old_R*(old_XYZvox-1) + old_T
      %
      %  We can use modified new_M1 & old_M1 to substitute new_M & old_M
      %      new_M1 * new_XYZvox       ==       old_M1 * old_XYZvox
      %
      %  where: M1 = M;   M1(:,4) = M(:,4) - sum(M(:,1:3),2);
      %  and:             M(:,4) == [T; 1] == sum(M1,2)
      %
      %  Therefore:   old_XYZvox = old_M1 \ new_M1 * new_XYZvox;
      %
      %  Since we are traverse Z axis, and   new_XYZvox   is replaced
      %  by   mask_Z * new_XYZvox, the above formula can be rewritten
      %  as:    old_XYZvox = old_M1 \ new_M1 * mask_Z * new_XYZvox;
      %
      %  i.e. we find the mapping from new_XYZvox to old_XYZvox:
      %  M = old_M1 \ new_M1 * mask_Z;
      %
      %  First, compute modified old_M1 & new_M1
      %
      old_M1 = old_M;   old_M1(:,4) = old_M(:,4) - sum(old_M(:,1:3),2);
      new_M1 = new_M;   new_M1(:,4) = new_M(:,4) - sum(new_M(:,1:3),2);

      %  Then, apply unit increment of mask_Z(3,4) to simulate the
      %  cursor movement
      %
      mask_Z(3,4) = z;

      %  Here is the mapping from new_XYZvox to old_XYZvox
      %
      M = old_M1 \ new_M1 * mask_Z;


      %switch method
      %case 1
      %   [new_img(:,:,z),origin(:,:,z,:)] = trilinear(old_img, new_dim, old_dim, M, bg);
      %case 2
         [new_img(:,:,z),origin(:,:,z,:)] = nearest_neighbor(old_img, new_dim, old_dim, M, bg);

      %end

   end;			% for z


   if ndims(old_img) == 2
      new_M(3,:) = [];
      new_M(:,3) = [];
   end

   return;					% affine



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [img_slice,origin_slice] = nearest_neighbor(img, dim1, dim2, M, bg)

   img_slice = zeros(dim1(1:2));
   origin_slice=zeros([dim1(1:2),3]);

   %  Dimension of transformed 3D volume
   %
   xdim1 = dim1(1);
   ydim1 = dim1(2);

   %  Dimension of original 3D volume
   %
   xdim2 = dim2(1);
   ydim2 = dim2(2);
   zdim2 = dim2(3);

   %  initialize new_Y accumulation
   %
   Y2X = 0;
   Y2Y = 0;
   Y2Z = 0;

   for y = 1:ydim1

      %  increment of new_Y accumulation
      %
      Y2X = Y2X + M(1,2);		% new_Y to old_X
      Y2Y = Y2Y + M(2,2);		% new_Y to old_Y
      Y2Z = Y2Z + M(3,2);		% new_Y to old_Z

      %  backproject new_Y accumulation and translation to old_XYZ
      %
      old_X = Y2X + M(1,4);
      old_Y = Y2Y + M(2,4);
      old_Z = Y2Z + M(3,4);

      for x = 1:xdim1

         %  accumulate the increment of new_X and apply it
         %  to the backprojected old_XYZ
         %
         old_X = M(1,1) + old_X  ;
         old_Y = M(2,1) + old_Y  ;
         old_Z = M(3,1) + old_Z  ;

         xi = round(old_X);
         yi = round(old_Y);
         zi = round(old_Z);

         %  within boundary of original image
         %
         if (	xi >= 1 & xi <= xdim2 & ...
		yi >= 1 & yi <= ydim2 & ...
		zi >= 1 & zi <= zdim2	)

            img_slice(x,y) = img(xi,yi,zi);
	    origin_slice(x,y,:)=[xi,yi,zi]';


         else
            img_slice(x,y) = bg;
	    origin_slice(x,y,:)=[0 0 0]';

         end	% if boundary

      end	% for x
   end		% for y

   return;					% nearest_neighbor


%--------------------------------------------------------------------

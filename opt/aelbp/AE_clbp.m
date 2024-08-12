% 
function [AECLBP_S,AECLBP_M,AECLBP_C] = AE_clbp(varargin)
% Check number of input arguments.
error(nargchk(1,5,nargin));

image=varargin{1};
d_image=double(image);

if nargin==1
    spoints=[-1 -1; -1 0; -1 1; 0 -1; -0 1; 1 -1; 1 0; 1 1];
    neighbors=8;
    mapping=0;
    mode='h';
end

if (nargin == 2) && (length(varargin{2}) == 1)
    error('Input arguments');
end

if (nargin > 2) && (length(varargin{2}) == 1)
    radius=varargin{2};
    neighbors=varargin{3};
    
    spoints=zeros(neighbors,2);

    % Angle step.
    a = 2*pi/neighbors;
    
    for i = 1:neighbors
        spoints(i,1) = -radius*sin((i-1)*a);
        spoints(i,2) = radius*cos((i-1)*a);
    end
    
    if(nargin >= 4)
        mapping=varargin{4};
        if(isstruct(mapping) && mapping.samples ~= neighbors)
            error('Incompatible mapping');
        end
    else
        mapping=0;
    end
    
    if(nargin >= 5)
        mode=varargin{5};
    else
        mode='h';
    end
end

if (nargin > 1) && (length(varargin{2}) > 1)
    spoints=varargin{2};
    neighbors=size(spoints,1);
    
    if(nargin >= 3)
        mapping=varargin{3};
        if(isstruct(mapping) && mapping.samples ~= neighbors)
            error('Incompatible mapping');
        end
    else
        mapping=0;
    end
    
    if(nargin >= 4)
        mode=varargin{4};
    else
        mode='h';
    end   
end

% Determine the dimensions of the input image.
[ysize xsize] = size(image);

miny=min(spoints(:,1));
maxy=max(spoints(:,1));
minx=min(spoints(:,2));
maxx=max(spoints(:,2));

% Block size, each LBP code is computed within a block of size bsizey*bsizex
bsizey=ceil(max(maxy,0))-floor(min(miny,0))+3;
bsizex=ceil(max(maxx,0))-floor(min(minx,0))+3;

% Coordinates of origin (0,0) in the block
origy=2-floor(min(miny,0));
origx=2-floor(min(minx,0));

% Minimum allowed size for the input image depends
% on the radius of the used operator.
if(xsize < bsizex || ysize < bsizey)
  error('Too small input image. Should be at least (2*radius+1) x (2*radius+1)');
end

% Calculate dx and dy;
dx = xsize - bsizex;
dy = ysize - bsizey;

% Fill the center pixel matrix C.
C = image(origy:origy+dy,origx:origx+dx);
d_C = double(C);

bins = 2^neighbors;

% Initialize the result matrix with zeros.
   AECLBP_M=zeros(dy+1,dx+1);
   AECLBP_C=zeros(dy+1,dx+1);
   AECLBP_S=zeros(dy+1,dx+1);
   NE_D=zeros(dy+1,dx+1);
%Compute the LBP code image

for i = 1:neighbors
  y = spoints(i,1)+origy;
  x = spoints(i,2)+origx;
  % Calculate floors, ceils and rounds for the x and y.
  ry = round(y);
  rx = round(x);
    NE1 = d_image(ry-1:ry-1+dy,rx-1:rx-1+dx);
    NE2 = d_image(ry-1:ry-1+dy,rx:rx+dx);
    NE3 = d_image(ry-1:ry-1+dy,rx+1:rx+1+dx);
    NE4 = d_image(ry:ry+dy,rx-1:rx-1+dx);
    NE5 = d_image(ry:ry+dy,rx+1:rx+1+dx);
    NE6 = d_image(ry+1:ry+1+dy,rx-1:rx-1+dx);
    NE7 = d_image(ry+1:ry+1+dy,rx:rx+dx);
    NE8 = d_image(ry+1:ry+1+dy,rx+1:rx+1+dx);
    NE_average = (NE1 + NE2 + NE3 + NE4 + NE5 + NE6 + NE7 + NE8)/8;
    A = NE_average >= d_C;
    Diff{i} = abs(NE_average-d_C);    
    MeanDiff(i) = mean(mean(Diff{i}));
    v = 2^(i-1);
    AECLBP_S = AECLBP_S +  v*A;
end
% Difference threshold for AECLBP_M
DiffThreshold = mean(MeanDiff);
% % compute AECLBP_M
for i=1:neighbors
  % Update the result matrix.
  v = 2^(i-1);
  AECLBP_M = AECLBP_M + v*(Diff{i}>=DiffThreshold);
end
% AECLBP_C
AECLBP_C = d_C>=mean(d_image(:));

%Apply mapping if it is defined
if isstruct(mapping)
    bins = mapping.num;
    for i = 1:size(AECLBP_S,1)
        for j = 1:size(AECLBP_S,2)
            AECLBP_M(i,j) = mapping.table(AECLBP_M(i,j)+1);
            AECLBP_S(i,j) = mapping.table(AECLBP_S(i,j)+1);
        end
    end
end

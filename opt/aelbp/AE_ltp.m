% 
function result = AE_ltp(varargin) % image,radius,neighbors,mapping,mode)
% Version for LTP

% Check number of input arguments.
error(nargchk(1,5,nargin));

image=varargin{1};
d_image=double(image);
%%%t thershold
t=5;

if nargin==1
    spoints=[-1 -1; -1 0; -1 1; 0 -1; -0 1; 1 -1; 1 0; 1 1];
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
%     spoints(i,1) = -radius*sin((i-1)*a+pi*3/4);
%     spoints(i,2) = radius*cos((i-1)*a+pi*3/4);
    if(nargin >= 4)
        mapping=varargin{4};
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
    
    if(nargin >= 3)
        mapping=varargin{3};
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

neighbors=size(spoints,1);

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
% on the radius of the used LBP operator.
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
result1=zeros(dy+1,dx+1);
result2=zeros(dy+1,dx+1);
%Compute the LBP code image

for i = 1:neighbors
  y = spoints(i,1)+origy;
  x = spoints(i,2)+origx;
  % Calculate floors, ceils and rounds for the x and y.
  fy = floor(y); cy = ceil(y); ry = round(y);
  fx = floor(x); cx = ceil(x); rx = round(x);
  % Check if interpolation is needed.
%   if (abs(x - rx) < 1e-6) && (abs(y - ry) < 1e-6)
    % Interpolation is not needed, use original datatypes
%     N = image(ry:ry+dy,rx:rx+dx);
    NE1 = d_image(ry-1:ry-1+dy,rx-1:rx-1+dx);
    NE2 = d_image(ry-1:ry-1+dy,rx:rx+dx);
    NE3 = d_image(ry-1:ry-1+dy,rx+1:rx+1+dx);
    NE4 = d_image(ry:ry+dy,rx-1:rx-1+dx);
    NE5 = d_image(ry:ry+dy,rx+1:rx+1+dx);
    NE6 = d_image(ry+1:ry+1+dy,rx-1:rx-1+dx);
    NE7 = d_image(ry+1:ry+1+dy,rx:rx+dx);
    NE8 = d_image(ry+1:ry+1+dy,rx+1:rx+1+dx);
    NE_average = (NE1 + NE2 + NE3 + NE4 + NE5 + NE6 + NE7 + NE8)/8;
%     A = NE_average >= d_C;
    D1 = NE_average >= (C+t);
    D2 = NE_average <= (C-t);
%   else
%     % Interpolation needed, use double type images 
%     ty = y - fy;
%     tx = x - fx;
% 
%     % Calculate the interpolation weights.
%     w1 = (1 - tx) * (1 - ty);
%     w2 =      tx  * (1 - ty);
%     w3 = (1 - tx) *      ty ;
%     w4 =      tx  *      ty ;
%     % Compute interpolated pixel values
%     N = w1*d_image(fy:fy+dy,fx:fx+dx) + w2*d_image(fy:fy+dy,cx:cx+dx) + ...
%         w3*d_image(cy:cy+dy,fx:fx+dx) + w4*d_image(cy:cy+dy,cx:cx+dx);
%     D1 = N >= (C+t);
%     D2 = N <= (C-t);
%   end  
  % Update the result matrix.
  v = 2^(i-1);
  result1 = result1 + v*D1;
  result2 = result2 + v*D2;
end

%Apply mapping if it is defined
if length(mapping) > 1
    bins = max(max(mapping)) + 1;
    for i = 1:size(result1,1)
        for j = 1:size(result1,2)
            result1(i,j) = mapping(result1(i,j)+1);
            result2(i,j) = mapping(result2(i,j)+1);
        end
    end
end

if (strcmp(mode,'h') || strcmp(mode,'hist') || strcmp(mode,'nh'))
    % Return with LBP histogram if mode equals 'hist'.
    result1=hist(result1(:),0:(bins-1));
    result2=hist(result2(:),0:(bins-1));
    result=[result1,result2];
end






function [resultp, resultn] = ltp2(varargin) % image,radius,neighbors,mapping,mode

% Check number of input arguments.
narginchk(1,6);
threshold=varargin{6};
image=varargin{1};
d_image=single(image);

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
[ysize, xsize] = size(image);



miny=min(spoints(:,1));
maxy=max(spoints(:,1));
minx=min(spoints(:,2));
maxx=max(spoints(:,2));

% Block size, each LBP code is computed within a block of size bsizey*bsizex
bsizey=ceil(max(maxy,0))-floor(min(miny,0))+1;
bsizex=ceil(max(maxx,0))-floor(min(minx,0))+1;

% Coordinates of origin (0,0) in the block
origy=1-floor(min(miny,0));
origx=1-floor(min(minx,0));

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
d_C = single(C);

bins = 2^neighbors;

% Initialize the result matrix with zeros.
resultp=zeros(dy+1,dx+1);
resultn=zeros(dy+1,dx+1);
%Compute the LBP code image

for i = 1:neighbors
  y = spoints(i,1)+origy;
  x = spoints(i,2)+origx;
  % Calculate floors, ceils and rounds for the x and y.
  fy = floor(y); cy = ceil(y); ry = round(y);
  fx = floor(x); cx = ceil(x); rx = round(x);
  % Check if interpolation is needed.
  if (abs(x - rx) < 1e-6) && (abs(y - ry) < 1e-6)
    % Interpolation is not needed, use original datatypes
    N = image(ry:ry+dy,rx:rx+dx);
    Dp = (N-C)>=threshold;
    Dn = (N-C)<= (-threshold);
  else
    % Interpolation needed, use double type images
    ty = y - fy;
    tx = x - fx;

    % Calculate the interpolation weights.
%     w1 = (1 - tx) * (1 - ty);
%     w2 =      tx  * (1 - ty);
%     w3 = (1 - tx) *      ty ;
%     w4 =      tx  *      ty ;
    % Compute interpolated pixel values
    a1=d_image(fy:fy+dy,fx:fx+dx);
    a2=d_image(fy:fy+dy,cx:cx+dx);
    a3=d_image(cy:cy+dy,fx:fx+dx);
    a4=d_image(cy:cy+dy,cx:cx+dx);
%     N = w1*d_image(fy:fy+dy,fx:fx+dx) + w2*d_image(fy:fy+dy,cx:cx+dx) + ...
%         w3*d_image(cy:cy+dy,fx:fx+dx) + w4*d_image(cy:cy+dy,cx:cx+dx);
    horInterp1 =a1*(1-tx)+a2*tx;
    horInterp2 =a3*(1-tx)+a4*tx;
    N=horInterp1*(1-ty)+horInterp2*ty;
    Dp = (N-d_C)>= threshold;
    Dn = (N-d_C)<=(-threshold);
  end 
  % Update the result matrix.
  v = 2^(i-1);
  resultp = resultp + v*Dp;
  resultn = resultn + v*Dn;
end

%Apply mapping if it is defined
if isstruct(mapping)
    bins = mapping.num;
    for i = 1:size(resultp,1)
        for j = 1:size(resultp,2)
            resultp(i,j) = mapping.table(resultp(i,j)+1);
            resultn(i,j) = mapping.table(resultn(i,j)+1);
        end
    end
end

if (strcmp(mode,'h') || strcmp(mode,'hist') || strcmp(mode,'nh'))
    % Return with LBP histogram if mode equals 'hist'.
    resultp=hist(resultp(:),0:(bins-1));
    resultn=hist(resultn(:),0:(bins-1));
%     resultpnh=histogram2(resultp(:),resultn(:),0:(bins-1),0:(bins-1));
%     resultpn=resultpnh.Values;
    result=[resultp resultn];
    if (strcmp(mode,'nh'))
        result=result/sum(result);
    end
end

end

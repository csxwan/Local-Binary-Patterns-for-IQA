function [ result ] = cldp( radius,n_sampling,img )
I=img;
mapping=getmapping(n_sampling,'riu2');
result=clderivativep(I,radius,n_sampling,mapping,'h');

end



function result = clderivativep(varargin) % image,radius,neighbors,mapping,mode)
% Version 0.3.2
% Authors: Marko Heikkil?and Timo Ahonen

% Changelog
% Version 0.3.2: A bug fix to enable using mappings together with a
% predefined spoints array
% Version 0.3.1: Changed MAPPING input to be a struct containing the mapping
% table and the number of bins to make the function run faster with high number
% of sampling points. Lauge Sorensen is acknowledged for spotting this problem.


% Check number of input arguments.
narginchk(1,5);

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
    if radius<=1
        error('radius too small');
    end
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
resultS=zeros(dy+1,dx+1);
resultD=zeros(dy+1,dx+1);
resultM=zeros(dy+1,dx+1);
c=zeros(dy+1,dx+1);

%Calculate average
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
  else
    % Interpolation needed, use double type images
    ty = y - fy;
    tx = x - fx;

%     % Calculate the interpolation weights.
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

  end 
    N=double(N);
    D = abs(N-d_C);
    c=c+D;
end
mc=mean(mean(c/neighbors));
mg=mean(mean(d_C));

%Compute the CLDP code image
for i = 1:neighbors
  y = spoints(i,1)+origy;
  x = spoints(i,2)+origx;
  y0 = spoints(i,1)*(radius-1)/radius+origy;
  x0 = spoints(i,2)*(radius-1)/radius+origx;
  % Calculate floors, ceils and rounds for the x and y.
  fy = floor(y); cy = ceil(y); ry = round(y);
  fx = floor(x); cx = ceil(x); rx = round(x);
  fy0 = floor(y0); cy0 = ceil(y0); ry0 = round(y0);
  fx0 = floor(x0); cx0 = ceil(x0); rx0 = round(x0);
  % Check if interpolation is needed.
  if (abs(x - rx) < 1e-6) && (abs(y - ry) < 1e-6)
    % Interpolation is not needed, use original datatypes
    N = image(ry:ry+dy,rx:rx+dx);
    N0 = image(ry0:ry0+dy,rx0:rx0+dx);
  else
    % Interpolation needed, use double type images
    ty = y - fy;
    tx = x - fx;
    ty0 = y0 - fy0;
    tx0 = x0 - fx0;
%     % Calculate the interpolation weights.
%     w1 = (1 - tx) * (1 - ty);
%     w2 =      tx  * (1 - ty);
%     w3 = (1 - tx) *      ty ;
%     w4 =      tx  *      ty ;
    % Compute interpolated pixel values
    a1=d_image(fy:fy+dy,fx:fx+dx);
    a2=d_image(fy:fy+dy,cx:cx+dx);
    a3=d_image(cy:cy+dy,fx:fx+dx);
    a4=d_image(cy:cy+dy,cx:cx+dx);
    a10=d_image(fy0:fy0+dy,fx0:fx0+dx);
    a20=d_image(fy0:fy0+dy,cx0:cx0+dx);
    a30=d_image(cy0:cy0+dy,fx0:fx0+dx);
    a40=d_image(cy0:cy0+dy,cx0:cx0+dx);
%     N = w1*d_image(fy:fy+dy,fx:fx+dx) + w2*d_image(fy:fy+dy,cx:cx+dx) + ...
%         w3*d_image(cy:cy+dy,fx:fx+dx) + w4*d_image(cy:cy+dy,cx:cx+dx);
    horInterp1 =a1*(1-tx)+a2*tx;
    horInterp2 =a3*(1-tx)+a4*tx;
    horInterp10 =a10*(1-tx0)+a20*tx0;
    horInterp20 =a30*(1-tx0)+a40*tx0;
    N=horInterp1*(1-ty)+horInterp2*ty;
    N0=horInterp10*(1-ty0)+horInterp20*ty0;
  end 
  N=double(N);
  S = N>=d_C;
  S0 = N0>=d_C;
  D=xor(S,S0);
  M = abs(N-d_C)>=mc;
  % Update the result matrix.
  v = 2^(i-1);
  resultS = resultS + v*S;
  resultD = resultD + v*D;
  resultM = resultM + v*M;
end
resultC=d_C>=mg;
%Apply mapping if it is defined
if isstruct(mapping)
    bins = mapping.num;
    for i = 1:size(resultS,1)
        for j = 1:size(resultS,2)
            resultS(i,j) = mapping.table(resultS(i,j)+1);
            resultD(i,j) = mapping.table(resultD(i,j)+1);
            resultM(i,j) = mapping.table(resultM(i,j)+1);
        end
    end
end

if (strcmp(mode,'h') || strcmp(mode,'hist') || strcmp(mode,'nh'))
    % Return with LBP histogram if mode equals 'hist'.
    resultS=hist(resultS(:),0:(bins-1));
    resultD=hist(resultD(:),0:(bins-1));
    hresultM_C=histogram2(resultC(:),resultM(:),0:2,0:(bins));
    resultM_C=hresultM_C.Values;
    result=[resultS resultD resultM_C(:)'];
    if (strcmp(mode,'nh'))
        result=result/sum(result);
    end
end

end

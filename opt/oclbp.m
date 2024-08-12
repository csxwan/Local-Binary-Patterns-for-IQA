function [result ] = oclbp( radius,n_sampling,img  )
I=img;
mapping=0;
result=oclbp2(I,radius,n_sampling,mapping,'h');


end


function result = oclbp2(varargin) % image,radius,neighbors,mapping,mode)
% Version 0.3.2
% Authors: Marko Heikkil?and Timo Ahonen

% Check number of input arguments.
narginchk(1,5);

image=varargin{1};
imageR=image(:,:,1);
imageG=image(:,:,2);
imageB=image(:,:,3);
d_image=single(image);
d_imageR=single(image(:,:,1));
d_imageG=single(image(:,:,2));
d_imageB=single(image(:,:,3));
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
[ysize, xsize,~] = size(image);



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
CR = imageR(origy:origy+dy,origx:origx+dx);
CG = imageG(origy:origy+dy,origx:origx+dx);
CB = imageB(origy:origy+dy,origx:origx+dx);
d_C = single(C);
d_CR = single(CR);
d_CG = single(CG);
d_CB = single(CB);

bins = 2^neighbors;

% Initialize the result matrix with zeros.
resultR=zeros(dy+1,dx+1);
resultG=zeros(dy+1,dx+1);
resultB=zeros(dy+1,dx+1);
resultRG=zeros(dy+1,dx+1);
resultRB=zeros(dy+1,dx+1);
resultGB=zeros(dy+1,dx+1);
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
    NR = imageR(ry:ry+dy,rx:rx+dx);
    NG = imageG(ry:ry+dy,rx:rx+dx);
    NB = imageB(ry:ry+dy,rx:rx+dx);
    DR = NR >=CR;
    DG = NG >=CG;
    DB = NB >=CB;
    DRG = NG>=CR;
    DRB = NB>=CR;
    DGB = NB>=CG;
  else
    % Interpolation needed, use double type images
    ty = y - fy;
    tx = x - fx;
    % Compute interpolated pixel values
    a1=d_image(fy:fy+dy,fx:fx+dx);
    a2=d_image(fy:fy+dy,cx:cx+dx);
    a3=d_image(cy:cy+dy,fx:fx+dx);
    a4=d_image(cy:cy+dy,cx:cx+dx);
    
    a1R=d_imageR(fy:fy+dy,fx:fx+dx);
    a2R=d_imageR(fy:fy+dy,cx:cx+dx);
    a3R=d_imageR(cy:cy+dy,fx:fx+dx);
    a4R=d_imageR(cy:cy+dy,cx:cx+dx);
    
    a1G=d_imageG(fy:fy+dy,fx:fx+dx);
    a2G=d_imageG(fy:fy+dy,cx:cx+dx);
    a3G=d_imageG(cy:cy+dy,fx:fx+dx);
    a4G=d_imageG(cy:cy+dy,cx:cx+dx);
    
    a1B=d_imageB(fy:fy+dy,fx:fx+dx);
    a2B=d_imageB(fy:fy+dy,cx:cx+dx);
    a3B=d_imageB(cy:cy+dy,fx:fx+dx);
    a4B=d_imageB(cy:cy+dy,cx:cx+dx);
%     N = w1*d_image(fy:fy+dy,fx:fx+dx) + w2*d_image(fy:fy+dy,cx:cx+dx) + ...
%         w3*d_image(cy:cy+dy,fx:fx+dx) + w4*d_image(cy:cy+dy,cx:cx+dx);
    horInterp1 =a1*(1-tx)+a2*tx;
    horInterp2 =a3*(1-tx)+a4*tx;
    
    horInterp1R =a1R*(1-tx)+a2R*tx;
    horInterp2R =a3R*(1-tx)+a4R*tx;
    
    horInterp1G =a1G*(1-tx)+a2G*tx;
    horInterp2G =a3G*(1-tx)+a4G*tx;
    
    horInterp1B =a1B*(1-tx)+a2B*tx;
    horInterp2B =a3B*(1-tx)+a4B*tx;
    
    N=horInterp1*(1-ty)+horInterp2*ty;
    NR=horInterp1R*(1-ty)+horInterp2R*ty;       
    NG=horInterp1G*(1-ty)+horInterp2G*ty;     
    NB=horInterp1B*(1-ty)+horInterp2B*ty;
    D = N>=d_C;
    DR = NR>=d_CR;
    DG = NG>=d_CG;
    DB = NB>=d_CB;
    DRG = NG>=d_CR;
    DRB = NB>=d_CR;
    DGB = NB>=d_CG;
  end 
  % Update the result matrix.
  v = 2^(i-1);
  resultR = resultR + v*DR;
  resultG = resultG + v*DG;
  resultB = resultB + v*DB;
  resultRG = resultRG + v*DRG;
  resultRB = resultRB + v*DRB;
  resultGB = resultGB + v*DGB;
end

%Apply mapping if it is defined
if isstruct(mapping)
    bins = mapping.num;
    for i = 1:size(resultR,1)
        for j = 1:size(resultR,2)
            resultR(i,j) = mapping.table(resultR(i,j)+1);
            resultG(i,j) = mapping.table(resultG(i,j)+1);
            resultB(i,j) = mapping.table(resultB(i,j)+1);
            resultRG(i,j) = mapping.table(resultRG(i,j)+1);
            resultRB(i,j) = mapping.table(resultRB(i,j)+1);
            resultGB(i,j) = mapping.table(resultGB(i,j)+1);
        end
    end
end

if (strcmp(mode,'h') || strcmp(mode,'hist') || strcmp(mode,'nh'))
    % Return with LBP histogram if mode equals 'hist'.
    resultR=hist(resultR(:),0:(bins-1));
    resultG=hist(resultG(:),0:(bins-1));
    resultB=hist(resultB(:),0:(bins-1));
    resultRG=hist(resultRG(:),0:(bins-1));
    resultRB=hist(resultRB(:),0:(bins-1));
    resultGB=hist(resultGB(:),0:(bins-1));
    result=[resultR resultG resultB resultRG resultRB resultGB];
    if (strcmp(mode,'nh'))
        result=result/sum(result);
    end
end

end

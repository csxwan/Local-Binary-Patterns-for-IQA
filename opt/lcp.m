function [ result ] = lcp( radius,n_sampling,img  )
I=img;
result=lcp2(I,radius,n_sampling);

end


function result = lcp2(varargin) % image,radius,neighbors

% Check number of input arguments.
narginchk(1,3);

image=varargin{1};
d_image=single(image);

if nargin==1
    spoints=[-1 -1; -1 0; -1 1; 0 -1; -0 1; 1 -1; 1 0; 1 1];
    neighbors=8;
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
end

if (nargin > 1) && (length(varargin{2}) > 1)
    spoints=varargin{2};
    neighbors=size(spoints,1);
end

% Determine the dimensions of the input image.
[ysize, xsize] = size(image);
I=image(radius+1:ysize-radius,radius+1:xsize-radius);


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
result1=zeros(dy+1,dx+1);
N=zeros(dy+1,dx+1,neighbors);

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
    N(:,:,i) = image(ry:ry+dy,rx:rx+dx);
    D = N(:,:,i) >=C;
  else
    % Interpolation needed, use double type images
    ty = y - fy;
    tx = x - fx;
    % Compute interpolated pixel values
    a1=d_image(fy:fy+dy,fx:fx+dx);
    a2=d_image(fy:fy+dy,cx:cx+dx);
    a3=d_image(cy:cy+dy,fx:fx+dx);
    a4=d_image(cy:cy+dy,cx:cx+dx);
%     N = w1*d_image(fy:fy+dy,fx:fx+dx) + w2*d_image(fy:fy+dy,cx:cx+dx) + ...
%         w3*d_image(cy:cy+dy,fx:fx+dx) + w4*d_image(cy:cy+dy,cx:cx+dx);
    horInterp1 =a1*(1-tx)+a2*tx;
    horInterp2 =a3*(1-tx)+a4*tx;
    N(:,:,i)=horInterp1*(1-ty)+horInterp2*ty;
    D = N(:,:,i)>=d_C;
  end 
  % Update the result matrix.
  v = 2^(i-1);
  result1 = result1 + v*D;
end
  resultlbp=result1;
  resulthist=hist(result1(:),0:(bins-1));
  result=[];
  [xx,yy]=size(resultlbp);
for n=1:bins
    temp=1;
    c=zeros(resulthist(n),1);
    V=zeros(resulthist(n),neighbors);
for i=1:xx
    for j=1:yy
      if resultlbp(i,j)==n-1
          c(temp,1)=I(i,j);
          V(temp,:)=N(i,j,:);
          temp=temp+1;
      end
    end
end
A=(V'*V)\(V'*c);
H1=fft(A);
H=abs(H1);
result=[result H' resulthist(n)];
end

end


function [ result ] = lgp( rad,n_sampling,img )
img=double(img);
N=1;
sigma1=0.5;
sigma2=2.5;

guass1x=guasskernelx(N,sigma1);
guass1y=guasskernely(N,sigma1);
guass2x=guasskernelx(N,sigma2);
guass2y=guasskernely(N,sigma2);

Gx1=imfilter(img, guass1x,'replicate');
Gx2=imfilter(img, guass2x,'replicate');
Gy1=imfilter(img, guass1y,'replicate');
Gy2=imfilter(img, guass2y,'replicate');

Gx1=double(Gx1);
Gy1=double(Gy1);
Gx2=double(Gx2);
Gy2=double(Gy2);

A1=sqrt(Gx1.^2+Gy1.^2);
A2=sqrt(Gx2.^2+Gy2.^2);
[m,n,h]=size(Gx1);
O1=zeros(m,n,h);
O2=zeros(m,n,h);

for H=1:h
for i=1:m
    for j=1:n
        if Gy1(i,j,H)~=0
           O1(i,j,H)=(atan(Gx1(i,j,H)/Gy1(i,j,H))+pi/2)/pi*360;
        end
        if Gy2(i,j,H)~=0
           O2(i,j,H)=(atan(Gx2(i,j,H)/Gy2(i,j,H))+pi/2)/pi*360;
        end
    end
end
end

map=getmapping(n_sampling,'riu2');
LGPA1=lbp2(A1,rad,n_sampling,map,'0');
LGPA2=lbp2(A2,rad,n_sampling,map,'0');

map=getmapping(n_sampling,'riu2');
LGPO1=lbpO(O1,rad,n_sampling,map,0);
LGPO2=lbpO(O2,rad,n_sampling,map,0);
 
hAO1=histogram2(LGPA1,LGPO1,0:n_sampling+2,0:n_sampling+2);
PAO1=hAO1.Values;
hAO2=histogram2(LGPA2,LGPO2,0:n_sampling+2,0:n_sampling+2);
PAO2=hAO2.Values;

resultA1=zeros(1,n_sampling+2);
resultO1=zeros(1,n_sampling+2);
resultA2=zeros(1,n_sampling+2);
resultO2=zeros(1,n_sampling+2);

 
for m=1:n_sampling+2
    for n=1:n_sampling+2
             resultA1(1,m)=resultA1(1,m)+PAO1(m,n)/sum(PAO1(:,n));
             resultA2(1,m)=resultA2(1,m)+PAO2(m,n)/sum(PAO2(:,n));
             resultO1(1,n)=resultO1(1,n)+PAO1(m,n)/sum(PAO1(m,:));
             resultO2(1,n)=resultO2(1,n)+PAO2(m,n)/sum(PAO2(m,:));

    end
end
result=[resultA1 resultO1 resultA2 resultO2]/(n_sampling+2);
end

function result = lbpO(varargin) % image,radius,neighbors,mapping,mode)
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
K=4;
C=getQ(C,K);
d_C = single(C);
d_C=getQ(d_C,K);

bins = 2^neighbors;

% Initialize the result matrix with zeros.
result=zeros(dy+1,dx+1);

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
    N=getQ(N,K);
    D = N == C; 
  else
    % Interpolation needed, use double type images 
    ty = y - fy;
    tx = x - fx;

    % Calculate the interpolation weights.
    w1 = (1 - tx) * (1 - ty);
    w2 =      tx  * (1 - ty);
    w3 = (1 - tx) *      ty ;
    w4 =      tx  *      ty ;
    % Compute interpolated pixel values
    N = w1*d_image(fy:fy+dy,fx:fx+dx) + w2*d_image(fy:fy+dy,cx:cx+dx) + ...
        w3*d_image(cy:cy+dy,fx:fx+dx) + w4*d_image(cy:cy+dy,cx:cx+dx);
   N=getQ(N,K);
    D = N == d_C; 
  end  
  % Update the result matrix.
  v = 2^(i-1);
  result = result + v*D;
end

%Apply mapping if it is defined
if isstruct(mapping)
    bins = mapping.num;
    for i = 1:size(result,1)
        for j = 1:size(result,2)
            result(i,j) = mapping.table(result(i,j)+1);
        end
    end
end

if (strcmp(mode,'h') || strcmp(mode,'hist') || strcmp(mode,'nh'))
    % Return with LBP histogram if mode equals 'hist'.
    result=hist(result(:),0:(bins-1));
    if (strcmp(mode,'nh'))
        result=result/sum(result);
    end
else
    %Otherwise return a matrix of unsigned integers
    if ((bins-1)<=intmax('uint8'))
        result=uint8(result);
    elseif ((bins-1)<=intmax('uint16'))
        result=uint16(result);
    else
        result=uint32(result);
    end
end

end

function Q=getQ(O,K)
Q=O;
[m,n]=size(O);
for q=1:K
 for i=1:m
     for j=1:n
         if O(i,j)>=360*(q-1)/K && O(i,j)<360*q/K 
              Q(i,j)=q;
         end
     end
 end
end

end
function mapping = getmapping(samples,mappingtype)
% Version 0.1.1
% Authors: Marko Heikkil?and Timo Ahonen

% Changelog
% 0.1.1 Changed output to be a structure
% Fixed a bug causing out of memory errors when generating rotation
% invariant mappings with high number of sampling points.
% Lauge Sorensen is acknowledged for spotting this problem.

 

table = 0:2^samples-1;
% index   = 0;

if strcmp(mappingtype,'riu2') 
  newMax = samples + 2;%Uniform & Rotation invariant
  for i = 0:2^samples - 1
    j = bitset(bitshift(i,1,'uint8'),1,bitget(i,samples)); %rotate left
    numt = sum(bitget(bitxor(i,j),1:samples));
    if numt <= 2
      table(i+1) = sum(bitget(i,1:samples));
    else
      table(i+1) = samples+1;
    end
  end
end
mapping.table=table;
mapping.samples=samples;
mapping.num=newMax;
end

function H=guasskernelx(N,sigma)
N_row = 2*N+1;
H = [];                                        %求高斯模板H
for i=1:N_row
    for j=1:N_row
        fenzi=double((i-N-1)^2+(j-N-1)^2);
        H(i,j)=(i-N-1)*exp(-fenzi/(2*sigma*sigma))/(-2*pi*sigma^4);
    end
end
S=sum(sum((H(1:N,:))));
H=H/S;
% H=H/sum(H(:));   
end

function H=guasskernely(N,sigma)
N_row = 2*N+1;
H = [];                                        %求高斯模板H
for i=1:N_row
    for j=1:N_row
        fenzi=double((i-N-1)^2+(j-N-1)^2);
        H(i,j)=(j-N-1)*exp(-fenzi/(2*sigma*sigma))/(-2*pi*sigma^4);
    end
end
S=sum(sum((H(:,1:N))));
H=H/S;
% H=H/sum(H(:));   
end
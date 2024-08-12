function [ result] = pi_lbp( radius,n_sampling,img )
I=img;
F2=[-1,1];
F3=[0,-1,1;
    -2,1,1;
    1,-2,1;
    1,1,-2];
F4=[0,0,-1,1;
    0,-2,1,1;
    0,1,-2,1;
    0,1,1,-2;
    1,-1,-1,1;
    1,-1,1,-1;
    1,1,-1,-1];

mapping=getmapping(n_sampling,'riu2');

result1=pilbp(I,radius,n_sampling,mapping,'h',F2,1);
temp=1;
for m=1:4
    for n=2:4
     result2(temp,:)=pilbp(I,radius,n_sampling,mapping,'h',F3(m,:),n);
     temp=temp+1;
    end
end
result2=result2(:)';
temp=1;
for m=1:7
    for n=5:7
     result3(temp,:)=pilbp(I,radius,n_sampling,mapping,'h',F4(m,:),n);
     temp=temp+1;
    end
end
result3=result3(:)';
result=[result1 result2 result3];
end

function result = pilbp(varargin) % image,radius,neighbors,mapping,mode)
% Version 0.3.2
% Authors: Marko Heikkil?and Timo Ahonen

% Changelog
% Version 0.3.2: A bug fix to enable using mappings together with a
% predefined spoints array
% Version 0.3.1: Changed MAPPING input to be a struct containing the mapping
% table and the number of bins to make the function run faster with high number
% of sampling points. Lauge Sorensen is acknowledged for spotting this problem.


% Check number of input arguments.

narginchk(1,7);

image=varargin{1};
d_image=single(image);
f=varargin{6};
g=varargin{7};
K=length(f);
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
    spoints1=zeros(neighbors,2);
    spoints2=zeros(neighbors,2);
    % Angle step.
    a = 2*pi/neighbors;
    a45=pi/4;
    for i = 1:neighbors
        spoints(i,1) = -radius*sin((i-1)*a);
        spoints(i,2) = radius*cos((i-1)*a);
        spoints1(i,1) = -radius*sin((i-1)*a-a45);
        spoints1(i,2) = radius*cos((i-1)*a-a45);
        spoints2(i,1) = -radius*sin((i-1)*a+a45);
        spoints2(i,2) = radius*cos((i-1)*a+a45);
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
bsizey=(K-1)*(ceil(max(maxy,0))-floor(min(miny,0))+1);
bsizex=(K-1)*(ceil(max(maxx,0))-floor(min(minx,0))+1);

% Coordinates of origin (0,0) in the block
origy=(K-1)*(1-floor(min(miny,0)));
origx=(K-1)*(1-floor(min(minx,0)));

% Minimum allowed size for the input image depends
% on the radius of the used LBP operator.
if(xsize < bsizex || ysize < bsizey)
  error('Too small input image. Should be at least (2*radius+1) x (2*radius+1)');
end

% Calculate dx and dy;
dx = xsize - bsizex;
dy = ysize - bsizey;

% Fill the center pixel matrix C.
% C = image(origy:origy+dy,origx:origx+dx);
% d_C = single(C);
D = zeros(dy+1,dx+1);
bins = 2^neighbors;

% Initialize the result matrix with zeros.
result=zeros(dy+1,dx+1);

%Compute the LBP code image
y=1:K;
x=1:K;

for i = 1:neighbors
    for k=1:K
     if k==1
          y(k) = origy;
          x(k) = origx;
     elseif k==K
         if g==2||g==5
          y(k) = spoints(i,1)+y(k-1);
          x(k)= spoints(i,2)+x(k-1); 
         elseif g==4||g==7
          y(k) = spoints1(i,1)+y(k-1);
          x(k)= spoints1(i,2)+x(k-1);  
         else
          y(k) = spoints2(i,1)+y(k-1);
          x(k)= spoints2(i,2)+x(k-1);  
         end
     else
          y(k) = spoints(i,1)+y(k-1);
          x(k)= spoints(i,2)+x(k-1);
     end   
          % Calculate floors, ceils and rounds for the x and y.
          fy = floor(y(k)); cy = ceil(y(k)); ry = round(y(k));
          fx = floor(x(k)); cx = ceil(x(k)); rx = round(x(k));
          % Check if interpolation is needed.
          if (abs(x(k) - rx) < 1e-6) && (abs(y(k) - ry) < 1e-6)
            % Interpolation is not needed, use original datatypes
            N = image(ry:ry+dy,rx:rx+dx);
            N=double(N);
            D=D+N*f(k);
          else
            % Interpolation needed, use double type images
            ty = y(k) - fy;
            tx = x(k) - fx;

%             % Calculate the interpolation weights.
%             w1 = (1 - tx) * (1 - ty);
%             w2 =      tx  * (1 - ty);
%             w3 = (1 - tx) *      ty ;
%             w4 =      tx  *      ty ;
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
            D=D+N*f(k);
          end 
     
    end
          D=D>=0;
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





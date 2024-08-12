function [ result ] = dlbp2(radius,n_sampling,img,size1)
% N=32;
% H1=circlegabor(N,2,0.414);
% H2=circlegabor(N,3.17,0.2612);
% H3=circlegabor(N,5.04,0.1643);
% H4=circlegabor(N,8,0.1035);
% % padSize = (N-1)/2;
% % image=rgb2gray(img);
% % a = padarray(image,padSize,'replicate');
% % A = fft2(a);
A=(img);
% I1=imfilter(A,H1,'same');
% % I11=I1(:);
% I2=imfilter(A,H2);
% I3=imfilter(A,H3);
% % I13=I3(:);
% I4=imfilter(A,H4);
A=rgb2gray(A);
A=imresize(A,size1);
A=double(A);
I1=imgaborf(A,2,0.414);
I2=imgaborf(A,3.17,0.2612);
I3=imgaborf(A,5.04,0.1643);
I4=imgaborf(A,8,0.1035);
I=I1+I2+I3+I4;

I1=abs(I1./I);
I2=abs(I2./I);
I3=abs(I3./I);
I4=abs(I4./I);
result2=[I1(:)' I2(:)' I3(:)' I4(:)'];
result1 = Dlbp(img,radius,n_sampling,'h');
result=[ result2 result1];
end



function result = Dlbp(varargin) % image,radius,neighbors,mapping,mode)
% Version 0.3.2
% Authors: Marko Heikkil?and Timo Ahonen

% Changelog
% Version 0.3.2: A bug fix to enable using mappings together with a
% predefined spoints array
% Version 0.3.1: Changed MAPPING input to be a struct containing the mapping
% table and the number of bins to make the function run faster with high number
% of sampling points. Lauge Sorensen is acknowledged for spotting this problem.


% Check number of input arguments.
narginchk(1,4);

image=varargin{1};
% k80=varargin{6};
d_image=single(image);

if nargin==1
    spoints=[-1 -1; -1 0; -1 1; 0 -1; -0 1; 1 -1; 1 0; 1 1];
    neighbors=8;
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
        mode=varargin{4};
    else
        mode='h';
    end
end

if (nargin > 1) && (length(varargin{2}) > 1)
    spoints=varargin{2};
    neighbors=size(spoints,1);
    
    
    if(nargin >= 3)
        mode=varargin{3};
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
    D = N >=C;
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
    a1=d_image(fy:fy+dy,fx:fx+dx);
    a2=d_image(fy:fy+dy,cx:cx+dx);
    a3=d_image(cy:cy+dy,fx:fx+dx);
    a4=d_image(cy:cy+dy,cx:cx+dx);
%     N = w1*d_image(fy:fy+dy,fx:fx+dx) + w2*d_image(fy:fy+dy,cx:cx+dx) + ...
%         w3*d_image(cy:cy+dy,fx:fx+dx) + w4*d_image(cy:cy+dy,cx:cx+dx);
    horInterp1 =a1*(1-tx)+a2*tx;
    horInterp2 =a3*(1-tx)+a4*tx;
    N=horInterp1*(1-ty)+horInterp2*ty;
    D = N>=d_C;
  end 
  % Update the result matrix.
  v = 2^(i-1);
  result = result + v*D;
end


if (strcmp(mode,'h') || strcmp(mode,'hist') || strcmp(mode,'nh'))
    % Return with LBP histogram if mode equals 'hist'.
    result=hist(result(:),0:(bins-1));
    result=sort(result,2,'descend');
%     result=result(1:k80);
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


function H=circlegabor(imageSize,F,sigma,classA)
    M = imageSize(1);
    N = imageSize(2);
    u = cast(createFrequencyVector(N),classA);
    v = cast(createFrequencyVector(M),classA);
    [U,V] = meshgrid(u,v);                                  
    fenzi=double((U).^2+(V).^2);
    H=exp(-2.*((fenzi).^0.5-F).^2.*(pi*sigma)^2)./(sqrt(2*pi)*sigma);  
 end

function ima=imgaborf(I,f,sigma)
    a=I;
    outSize = size(a);
    sizeLargestKernel=[3 3];
    padSize = (sizeLargestKernel-1)/2;
    a = padarray(a,padSize,'replicate');
    sizeAPadded = size(a);
    A = fft2(a);
    H=circlegabor(sizeAPadded,f,sigma,class(A));
    outPadded = ifft2(A .* ifftshift(H));
    out= unpadSlice(outPadded,padSize,outSize);
    ima=abs(out);
end

function u = createFrequencyVector(N)

if mod(N,2)
    u = linspace(-0.5+1/(2*N),0.5-1/(2*N),N);
else
    u = linspace(-0.5,0.5-1/N,N); 
end 

end

function outTrimmed = unpadSlice(out,padSize,outSize)

    start = padSize+1;
    stop = start+outSize-1;
    outTrimmed = out(start(1):stop(1),start(2):stop(2));

end
function [ result ] = mltp(k,n_sampling,img)
I=img;
T=1;
mapping=getmap(n_sampling,'riu3');
result=medianltp(I,k,n_sampling,mapping,'h',T);
end


function result = medianltp(varargin) % image,radius,neighbors,mapping,mode,threshold
% Version 0.3.2
% Authors: Marko Heikkil?and Timo Ahonen

% Changelog
% Version 0.3.2: A bug fix to enable using mappings together with a
% predefined spoints array
% Version 0.3.1: Changed MAPPING input to be a struct containing the mapping
% table and the number of bins to make the function run faster with high number
% of sampling points. Lauge Sorensen is acknowledged for spotting this problem.


% Check number of input arguments.
narginchk(1,6);

image=varargin{1};
d_image=single(image);
k=varargin{2};
T=varargin{6};
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
    neighbors=varargin{3};
    radius=(2*k+1)*neighbors/8;
    n=(2*k+1)*neighbors;
    spoints=zeros(n,2);

    % Angle step.
    a = 2*pi/n;
    
    for i = 1:n
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

bins = 3^neighbors;

% Initialize the result matrix with zeros.
resultR=zeros(dy+1,dx+1);
resultM=zeros(dy+1,dx+1);
% resultC=zeros(dy+1,dx+1);
N=zeros(dy+1,dx+1,2*k+1);
c=zeros(dy+1,dx+1);

%求全图的差值平均值
id=1;
for i = 1:n
  y = spoints(i,1)+origy;
  x = spoints(i,2)+origx;
  % Calculate floors, ceils and rounds for the x and y.
  fy = floor(y); cy = ceil(y); ry = round(y);
  fx = floor(x); cx = ceil(x); rx = round(x);
  % Check if interpolation is needed.
  if (abs(x - rx) < 1e-6) && (abs(y - ry) < 1e-6)
    % Interpolation is not needed, use original datatypes
    N(:,:,id) = image(ry:ry+dy,rx:rx+dx);
    id=id+1;
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
    N(:,:,id)=horInterp1*(1-ty)+horInterp2*ty;
    id=id+1;
  end 
  if mod(i,2*k+1)==0
      M=median(N(:,:,:),3);
      D=abs( M-d_C);
      c=c+D;
      id=1;
  end
end
mc=mean(mean(c/neighbors));
mg=mean(mean(d_C));
%Compute the LBP code image
id=1;
ide=1;
for i = 1:n
  y = spoints(i,1)+origy;
  x = spoints(i,2)+origx;
  % Calculate floors, ceils and rounds for the x and y.
  fy = floor(y); cy = ceil(y); ry = round(y);
  fx = floor(x); cx = ceil(x); rx = round(x);
  % Check if interpolation is needed.
  if (abs(x - rx) < 1e-6) && (abs(y - ry) < 1e-6)
    % Interpolation is not needed, use original datatypes
    N(:,:,id) = image(ry:ry+dy,rx:rx+dx);
    id=id+1;
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
    N(:,:,id)=horInterp1*(1-ty)+horInterp2*ty;
    id=id+1;
%     D = N>=d_C;
  end 
  if mod(i,2*k+1)==0
      M=median(N(:,:,:),3);
      R = M-d_C>=-T;
      R2=M-d_C>=T;
      R1=R-R2;
      mM =abs(M-d_C)-mc>=-T;
      mM2=abs(M-d_C)-mc>=T;
      mM1=mM-mM2;
      % Update the result matrix.
      v = 3^(ide-1);
      resultR = resultR + v*R1+v*2*R2;
      resultM = resultM + v*mM1+v*2*mM2;
      id=1;
      ide=ide+1;
  end
end
mC=d_C-mg>=-T;
mC2=d_C-mg>=T;
mC1=mC-mC2;
resultC=mC1+2*mC2;
%Apply mapping if it is defined
if isstruct(mapping)
    bins = mapping.num;
    for i = 1:size(resultR,1)
        for j = 1:size(resultR,2)
            resultR(i,j) = mapping.table(resultR(i,j)+1);
            resultM(i,j) = mapping.table(resultM(i,j)+1);
        end
    end
end

if (strcmp(mode,'h') || strcmp(mode,'hist') || strcmp(mode,'nh'))
    % Return with LBP histogram if mode equals 'hist'.
%     resultR=hist(resultR(:),0:(bins-1));
%     resultM=hist(resultM(:),0:(bins-1));
%     resultC=hist(resultC(:),0:2);
    resultR0=resultR;
    resultR0(resultC==1)=-1;
    resultR0(resultC==2)=-1;
    resultM0=resultM;
    resultM0(resultC==1)=-1;
    resultM0(resultC==2)=-1;
    H0=histogram2(resultM0, resultR0,0:(bins),0:(bins));
    hH0=H0.Values;
    
    resultR1=resultR;
    resultR1(resultC==0)=-1;
    resultR1(resultC==2)=-1;
    resultM1=resultM;
    resultM1(resultC==0)=-1;
    resultM1(resultC==2)=-1;
    H1=histogram2(resultM1, resultR1,0:(bins),0:(bins));
    hH1=H1.Values;
    
    resultR2=resultR;
    resultR2(resultC==1)=-1;
    resultR2(resultC==0)=-1;
    resultM2=resultM;
    resultM2(resultC==1)=-1;
    resultM2(resultC==0)=-1;
    H2=histogram2(resultM2, resultR2,0:(bins),0:(bins));
    hH2=H2.Values;
    
    result=[hH0(:)' hH1(:)' hH2(:)'];
    if (strcmp(mode,'nh'))
        result=result/sum(result);
    end
end

end



% GETMAPPING returns a structure containing a mapping table for LBP codes.
%  MAPPING = GETMAPPING(SAMPLES,MAPPINGTYPE) returns a
%  structure containing a mapping table for
%  LBP codes in a neighbourhood of SAMPLES sampling
%  points. Possible values for MAPPINGTYPE are
%       'riu3' for uniform rotation-invariant LBP.
%
%  Example:
%       I=imread('rice.tif');
%       MAPPING=getmapping(16,'riu3');
%       LBPHIST=lbp(I,2,16,MAPPING,'hist');
%  Now LBPHIST contains a rotation-invariant uniform LBP
%  histogram in a (16,2) neighbourhood.
%
function mapping = getmap(samples,mappingtype)
% Version 0.1.1
% Authors: Marko Heikkil?and Timo Ahonen

% Changelog
% 0.1.1 Changed output to be a structure
% Fixed a bug causing out of memory errors when generating rotation
% invariant mappings with high number of sampling points.
% Lauge Sorensen is acknowledged for spotting this problem.

 

table = 0:3^samples-1;
newMax  = 0; %number of patterns in the resulting LBP code
% index   = 0;

% if strcmp(mappingtype,'u3') %Uniform 3
%   newMax = samples*(samples-1) + 3;
%   for i = 0:3^samples-1
% %       if(i==196)
% %          l=1;
% %       end
%     a=bitshift(i,1,'uint8');
%     b=bitget(i,samples);
%     j = bitset(bitshift(i,1,'uint8'),1,bitget(i,samples)); %rotate left
%     c=bitget(bitxor(i,j),1:samples);
%     numt = sum(bitget(bitxor(i,j),1:samples)); %number of 1->0 and
%                                                %0->1 transitions
%                                                %in binary string
%                                                %x is equal to the
%                                                %number of 1-bits in
%                                                %XOR(x,Rotate left(x))
%     if numt <= 2
%       table(i+1) = index;
%       index = index + 1;
%     else
%       table(i+1) = newMax - 1;
%     end
%   end
% end

% if strcmp(mappingtype,'ri') %Rotation invariant
%   tmpMap = zeros(3^samples,1) - 1;
%   for i = 0:3^samples-1
%     rm = i;
%     r  = i;
%     for j = 1:samples-1
%       r = bitset(bitshift(r,1,'uint8'),1,bitget(r,samples)); %rotate
%                                                              %left
%       if r < rm
%         rm = r;
%       end
%     end
%     if tmpMap(rm+1) < 0
%       tmpMap(rm+1) = newMax;
%       newMax = newMax + 1;
%     end
%     table(i+1) = tmpMap(rm+1);
%   end
% end

if strcmp(mappingtype,'riu3') %Uniform & Rotation invariant
  newMax = 3*samples;
  for i = 0:3^samples - 1
    i1 = dec2base(i,3,samples);
    i1=STR2DOUBLE(i1(:))';
    j=(circshift(i1,[0,1]));%rotate right
    
    A=bitxor(i1,j);
    numt = sum(A~=0);
    if numt <= 3&&sum(i1==1)<=1
        min=i;
       for n=1:7
          j1=circshift(i1,[0,n]);
          j1=num2str(j1);
          j1=strrep(j1, ' ', '');
          j2=base2dec(j1,3);
          if j2<min
              min=j2;
          end
       end
       i2=dec2base(min,3,samples);
       i2=str2num(i2(:))';
       if sum(i2==1)<1
          table(i+1) = sum(i2==2);
       elseif find(i2==1)==samples&&sum(i2==0)~=0
          table(i+1) = samples+sum(i1==0);  
       else
          table(i+1)= 2*samples-1+sum(i1==2);
       end
    else
      table(i+1) = 3*samples-1;
    end
  end
end

mapping.table=table;
mapping.samples=samples;
mapping.num=newMax;
end
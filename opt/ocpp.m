function [ result  ] = ocpp( radius,neighborsxy,neighborsxz,neighborsyz,img )
I=rgb2hsv(img);
mapping=0;
result=ocpp2(I,radius,neighborsxy,neighborsxz,neighborsyz,mapping,'h');


end


function result = ocpp2(varargin) % image,radius,neighbors,mapping,mode)
% Version 0.3.2
% Authors: Marko Heikkil?and Timo Ahonen

% Check number of input arguments.
narginchk(1,7);

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
    neighborsxy=varargin{3};
    neighborsxz=varargin{4};
    neighborsyz=varargin{5};
    
    spointsxy=zeros(neighborsxy,2);
    spointsxz=zeros(neighborsxz,2);
    spointsyz=zeros(neighborsyz,2);
    % Angle step.
    axy = 2*pi/neighborsxy;
    axz = 2*pi/neighborsxz;
    ayz = 2*pi/neighborsyz;
    for i = 1:neighborsxy
        spointsxy(i,1) = -radius*sin((i-1)*axy);
        spointsxy(i,2) = radius*cos((i-1)*axy);
    end
    for j = 1:neighborsxz
        spointsxz(j,1) = -sin((j-1)*axz);
        spointsxz(j,2) = radius*cos((j-1)*axz);
    end
    for k = 1:neighborsyz
        spointsyz(k,1) = -sin((k-1)*ayz);
        spointsyz(k,2) = radius*cos((k-1)*ayz);
    end
    
    if(nargin >= 6)
        mapping=varargin{6};
        if(isstruct(mapping) && mapping.samples ~= neighbors)
            error('Incompatible mapping');
        end
    else
        mapping=0;
    end
    
    if(nargin >= 7)
        mode=varargin{7};
    else
        mode='h';
    end
end

% if (nargin > 1) && (length(varargin{2}) > 1)
%     spoints=varargin{2};
%     neighbors=size(spoints,1);
%     
%     if(nargin >= 3)
%         mapping=varargin{3};
%         if(isstruct(mapping) && mapping.samples ~= neighbors)
%             error('Incompatible mapping');
%         end
%     else
%         mapping=0;
%     end
%     
%     if(nargin >= 4)
%         mode=varargin{4};
%     else
%         mode='h';
%     end   
% end
mappingxy=getmapping(neighborsxy,'u2');
mappingxz=getmapping(neighborsxz,'u2');
mappingyz=getmapping(neighborsyz,'u2');
% Determine the dimensions of the input image.
[ysize, xsize,zsize] = size(image);



minyxy=min(spointsxy(:,1));
maxyxy=max(spointsxy(:,1));
minxxy=min(spointsxy(:,2));
maxxxy=max(spointsxy(:,2));

% Block size, each LBP code is computed within a block of size bsizey*bsizex
bsizeyxy=ceil(max(maxyxy,0))-floor(min(minyxy,0))+1;
bsizexxy=ceil(max(maxxxy,0))-floor(min(minxxy,0))+1;

% Coordinates of origin (0,0) in the block
origyxy=1-floor(min(minyxy,0));
origxxy=1-floor(min(minxxy,0));

% Minimum allowed size for the input image depends
% on the radius of the used LBP operator.
if(xsize < bsizexxy || ysize < bsizeyxy)
  error('Too small input image. Should be at least (2*radius+1) x (2*radius+1)');
end

% Calculate dx and dy;
dxxy = xsize - bsizexxy;
dyxy = ysize - bsizeyxy;

% Fill the center pixel matrix C.
Cxy = image(origyxy:origyxy+dyxy,origxxy:origxxy+dxxy,:);
d_Cxy = single(Cxy);





binsxy = 2^neighborsxy;
binsxz = 2^neighborsxz;
binsyz = 2^neighborsyz;

% Initialize the result matrix with zeros.
resultxy=zeros(dyxy+1,dxxy+1,zsize);



%Compute the LBP code image

for i = 1:neighborsxy
  yxy = spointsxy(i,1)+origyxy;
  xxy = spointsxy(i,2)+origxxy;
  % Calculate floors, ceils and rounds for the x and y.
  fyxy = floor(yxy); cyxy = ceil(yxy); ryxy = round(yxy);
  fxxy = floor(xxy); cxxy = ceil(xxy); rxxy = round(xxy);
  % Check if interpolation is needed.
  if (abs(xxy - rxxy) < 1e-6) && (abs(yxy - ryxy) < 1e-6)
    % Interpolation is not needed, use original datatypes
    Nxy = image(ryxy:ryxy+dyxy,rxxy:rxxy+dxxy,:);
    Dxy = Nxy >=Cxy;
  else
    % Interpolation needed, use double type images
    tyxy = yxy - fyxy;
    txxy = xxy - fxxy;
    % Compute interpolated pixel values
    a1xy=d_image(fyxy:fyxy+dyxy,fxxy:fxxy+dxxy,:);
    a2xy=d_image(fyxy:fyxy+dyxy,cxxy:cxxy+dxxy,:);
    a3xy=d_image(cyxy:cyxy+dyxy,fxxy:fxxy+dxxy,:);
    a4xy=d_image(cyxy:cyxy+dyxy,cxxy:cxxy+dxxy,:);
%     N = w1*d_image(fy:fy+dy,fx:fx+dx) + w2*d_image(fy:fy+dy,cx:cx+dx) + ...
%         w3*d_image(cy:cy+dy,fx:fx+dx) + w4*d_image(cy:cy+dy,cx:cx+dx);
    horInterp1xy =a1xy*(1-txxy)+a2xy*txxy;
    horInterp2xy =a3xy*(1-txxy)+a4xy*txxy;
    Nxy=horInterp1xy*(1-tyxy)+horInterp2xy*tyxy;
    Dxy = Nxy>=d_Cxy;
  end 
  % Update the result matrix.
  vxy = 2^(i-1);
  resultxy = resultxy + vxy*Dxy;
end


minzxz=min(spointsxz(:,1));
maxzxz=max(spointsxz(:,1));
minxxz=min(spointsxz(:,2));
maxxxz=max(spointsxz(:,2));

% Block size, each LBP code is computed within a block of size bsizey*bsizex
bsizezxz=ceil(max(maxzxz,0))-floor(min(minzxz,0))+1;
bsizexxz=ceil(max(maxxxz,0))-floor(min(minxxz,0))+1;

% Coordinates of origin (0,0) in the block
origzxz=1-floor(min(minzxz,0));
origxxz=1-floor(min(minxxz,0));

% Minimum allowed size for the input image depends
% on the radius of the used LBP operator.
if(ysize < bsizexxz || zsize < bsizezxz)
  error('Too small input image. Should be at least (2*radius+1) x (2*radius+1)');
end

% Calculate dx and dy;
dxxz = ysize - bsizexxz;
dzxz = zsize - bsizezxz;

Cxz = image(origxxz:origxxz+dxxz,:,origzxz:origzxz+dzxz);
d_Cxz = single(Cxz);

resultxz=zeros(dxxz+1,xsize,dzxz+1);

for j = 1:neighborsxz
  zxz = spointsxz(j,1)+origzxz;
  xxz = spointsxz(j,2)+origxxz;
  % Calculate floors, ceils and rounds for the x and y.
  fzxz = floor(zxz); czxz = ceil(zxz); rzxz = round(zxz);
  fxxz = floor(xxz); cxxz = ceil(xxz); rxxz = round(xxz);
  % Check if interpolation is needed.
  if (abs(xxz - rxxz) < 1e-6) && (abs(zxz - rzxz) < 1e-6)
    % Interpolation is not needed, use original datatypes
    Nxz = image(rxxz:rxxz+dxxz,:,rzxz:rzxz+dzxz);
    Dxz = Nxz >=Cxz;
  else
    % Interpolation needed, use double type images
    tzxz = zxz - fzxz;
    txxz = xxz - fxxz;


    % Compute interpolated pixel values
    a1xz=d_image(fxxz:fxxz+dxxz,:,fzxz:fzxz+dzxz);
    a2xz=d_image(cxxz:cxxz+dxxz,:,fzxz:fzxz+dzxz);
    a3xz=d_image(fxxz:fxxz+dxxz,:,czxz:czxz+dzxz);
    a4xz=d_image(cxxz:cxxz+dxxz,:,czxz:czxz+dzxz);
%     N = w1*d_image(fy:fy+dy,fx:fx+dx) + w2*d_image(fy:fy+dy,cx:cx+dx) + ...
%         w3*d_image(cy:cy+dy,fx:fx+dx) + w4*d_image(cy:cy+dy,cx:cx+dx);
    horInterp1xz =a1xz*(1-txxz)+a2xz*txxz;
    horInterp2xz =a3xz*(1-txxz)+a4xz*txxz;
    Nxz=horInterp1xz*(1-tzxz)+horInterp2xz*tzxz;
    Dxz = Nxz>=d_Cxz;
  end 
  % Update the result matrix.
  vxz = 2^(j-1);
  resultxz = resultxz + vxz*Dxz;
end



minzyz=min(spointsyz(:,1));
maxzyz=max(spointsyz(:,1));
minyyz=min(spointsyz(:,2));
maxyyz=max(spointsyz(:,2));

% Block size, each LBP code is computed within a block of size bsizey*bsizex
bsizezyz=ceil(max(maxzyz,0))-floor(min(minzyz,0))+1;
bsizeyyz=ceil(max(maxyyz,0))-floor(min(minyyz,0))+1;

% Coordinates of origin (0,0) in the block
origzyz=1-floor(min(minzyz,0));
origyyz=1-floor(min(minyyz,0));

% Minimum allowed size for the input image depends
% on the radius of the used LBP operator.
if(xsize < bsizeyyz || zsize < bsizezyz)
  error('Too small input image. Should be at least (2*radius+1) x (2*radius+1)');
end

% Calculate dx and dy;
dyyz = xsize - bsizeyyz;
dzyz = zsize - bsizezyz;

Cyz = image(:,origyyz:origyyz+dyyz,origzyz:origzyz+dzyz);
d_Cyz = single(Cyz);

resultyz=zeros(ysize,dyyz+1,dzyz+1);

for k = 1:neighborsyz
  zyz = spointsyz(k,1)+origzyz;
  yyz = spointsyz(k,2)+origyyz;
  % Calculate floors, ceils and rounds for the x and y.
  fzyz = floor(zyz); czyz = ceil(zyz); rzyz = round(zyz);
  fyyz = floor(yyz); cyyz = ceil(yyz); ryyz = round(yyz);
  % Check if interpolation is needed.
  if (abs(yyz - ryyz) < 1e-6) && (abs(zyz - rzyz) < 1e-6)
    % Interpolation is not needed, use original datatypes
    Nyz = image(:,ryyz:ryyz+dyyz,rzyz:rzyz+dzyz);
    Dyz = Nyz >=Cyz;
  else
    % Interpolation needed, use double type images
    tzyz = zyz - fzyz;
    tyyz = yyz - fyyz;

    % Compute interpolated pixel values
    a1yz=d_image(:,fyyz:fyyz+dyyz,fzyz:fzyz+dzyz);
    a2yz=d_image(:,cyyz:cyyz+dyyz,fzyz:fzyz+dzyz);
    a3yz=d_image(:,fyyz:fyyz+dyyz,czyz:czyz+dzyz);
    a4yz=d_image(:,cyyz:cyyz+dyyz,czyz:czyz+dzyz);
%     N = w1*d_image(fy:fy+dy,fx:fx+dx) + w2*d_image(fy:fy+dy,cx:cx+dx) + ...
%         w3*d_image(cy:cy+dy,fx:fx+dx) + w4*d_image(cy:cy+dy,cx:cx+dx);
    horInterp1yz =a1yz*(1-tyyz)+a2yz*tyyz;
    horInterp2yz =a3yz*(1-tyyz)+a4yz*tyyz;
    Nyz=horInterp1yz*(1-tzyz)+horInterp2yz*tzyz;
    Dyz = Nyz>=d_Cyz;
  end 
  % Update the result matrix.
  vyz = 2^(k-1);
  resultyz = resultyz + vyz*Dyz;
end

resultxy1=resultxy(:,:,1);
resultxy2=resultxy(:,:,2);
resultxy3=resultxy(:,:,3);
[resultsizex,resultsizey]=size(resultxy2);
resultxz=resultxz(:,radius+1:radius+resultsizey);
resultyz=resultyz(radius+1:radius+resultsizex,:);


%Apply mapping if it is defined
if isstruct(mapping)
    binsxy = mappingxy.num;
    binsyz = mappingyz.num;
    binsxz = mappingxz.num;
    for i = 1:size(resultxy2,1)
        for j = 1:size(resultxy2,2)
            resultxy1(i,j) = mappingxy.table(resultxy1(i,j)+1);
            resultxy2(i,j) = mappingxy.table(resultxy2(i,j)+1);
            resultxy3(i,j) = mappingxy.table(resultxy3(i,j)+1);
            resultyz(i,j) = mappingyz.table(resultyz(i,j)+1);
            resultxz(i,j) = mappingxz.table(resultxz(i,j)+1);
        end
    end
end

if (strcmp(mode,'h') || strcmp(mode,'hist') || strcmp(mode,'nh'))
    % Return with LBP histogram if mode equals 'hist'.
    resultxy1=hist(resultxy1(:),0:(binsxy-1));
%     resultxy2=hist(resultxy2(:),0:(binsxy-1));
    resultxy3=hist(resultxy3(:),0:(binsxy-1));
%     resultxz=hist(resultxz(:),0:(binsxz-1));
%     resultyz=hist(resultyz(:),0:(binsyz-1));
    resultxy1=hist(resultxy1(:));
%     resultxy2=hist(resultxy2(:));
    resultxy3=hist(resultxy3(:));
%     resultxz=hist(resultxz(:));
%     resultyz=hist(resultyz(:));

%     CLBP_SM = [resultxz(:),resultyz(:)];
%     Hist3DSM = hist3(CLBP_SM,[mappingxz.num,mappingyz.num]);
%     CLBP_SMH = reshape(Hist3DSM,1,numel(Hist3DSM));

%     A=resultxy(:);
%     B=CLBP_SMH(:);
% %     CLBP_SMC = [A,B];
% %     Hist3D = hist3(CLBP_SMC,'Nbins',[mappingxy.num,mappingxz.num*mappingyz.num]);
%     Hist3D = histogram2(A,B,[mappingxy.num,mappingxz.num*mappingyz.num]);
%     CLBP_SMCH = reshape(Hist3D,1,numel(Hist3D));
%     result=CLBP_SMCH;
    
%     result=[resultxy1 resultxy2 resultxy3 CLBP_SMH];

   for i=0:mappingxz.num
    idx=find(resultxz==i);
    resultyz(idx)=resultyz(idx)+(i+1)*mappingyz.num;
   end
   CLBP_SM = [resultxy2(:),resultyz(:)];
   Hist3DSM = hist3(CLBP_SM,[mappingxy.num,mappingxz.num*mappingyz.num]);
   CLBP_SMH = reshape(Hist3DSM,1,numel(Hist3DSM));
   result=[resultxy1 resultxy3 CLBP_SMH];

%     result1=histogram2(resultxy(:),resultxz(:),0:(binsxy-1),0:(binsxz-1));
%     result1=result1.Values;
%     result1=result1(:);
%     result=histogram2(result1(:),resultyz(:),0:(binsyz-1),0:(binsyz-1));
%     result=result.Values;
    if (strcmp(mode,'nh'))
        result=result/sum(result);
    end
end

end

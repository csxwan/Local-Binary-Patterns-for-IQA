function [ result ] = dlbp( radius,n_sampling,img  )
I=img;
mapping=getmapping(n_sampling,'riu2');
result=Decorrelatedlbp(I,radius,n_sampling,mapping,'h');
end


function result = Decorrelatedlbp(varargin) % image,radius,neighbors,mapping,mode
 
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
    K=radius;
    spoints=zeros(neighbors,2);
    N=3*neighbors;
    % Angle step.
    a = 2*pi/N;
    
    for i = 1:N
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
    N=size(spoints,1);
    
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
resultd1=zeros(dy+1,dx+1);
resultd2=zeros(dy+1,dx+1);
resultd3=zeros(dy+1,dx+1);
resultM1=zeros(dy+1,dx+1);
resultM2=zeros(dy+1,dx+1);
resultM3=zeros(dy+1,dx+1);
% Initialize the required matrix with zeros.
N=zeros(dy+1,dx+1,K);
a1=zeros(dy+1,dx+1,K);
a1n=zeros(dy+1,dx+1,K);
a1p=zeros(dy+1,dx+1,K);
a2=zeros(dy+1,dx+1,K);
a2n=zeros(dy+1,dx+1,K);
a2p=zeros(dy+1,dx+1,K);
a3=zeros(dy+1,dx+1,K);
a3n=zeros(dy+1,dx+1,K);
a3p=zeros(dy+1,dx+1,K);
a4=zeros(dy+1,dx+1,K);
a4n=zeros(dy+1,dx+1,K);
a4p=zeros(dy+1,dx+1,K);
horInterp1=zeros(dy+1,dx+1,K);
horInterp2=zeros(dy+1,dx+1,K);
horInterp1n=zeros(dy+1,dx+1,K);
horInterp2n=zeros(dy+1,dx+1,K);
horInterp1p=zeros(dy+1,dx+1,K);
horInterp2p=zeros(dy+1,dx+1,K);
A1=zeros(dy+1,dx+1);
A2=zeros(dy+1,dx+1);
A3=zeros(dy+1,dx+1);
MM1=zeros(dy+1,dx+1);
MM2=zeros(dy+1,dx+1);
MM3=zeros(dy+1,dx+1);
% Initialize the required variables.
y=1:K;
x=1:K;
yn=1:K;
xn=1:K;
yp=1:K;
xp=1:K;
fy=1:K;
fx=1:K;
fyn=1:K;
fxn=1:K;
fyp=1:K;
fxp=1:K;
cy=1:K;
cx=1:K;
cyn=1:K;
cxn=1:K;
cyp=1:K;
cxp=1:K;
ry=1:K;
rx=1:K;
ryn=1:K;
rxn=1:K;
ryp=1:K;
rxp=1:K;
ty=1:K;
tx=1:K;
tyn=1:K;
txn=1:K;
typ=1:K;
txp=1:K;

%Compute the threshold
for i = 1:neighbors
  for k=1:radius  
      y(k) = spoints(i,1)*(k/radius)+origy;
      x(k) = spoints(i,2)*(k/radius)+origx;
      %Dense sampling
      if i==neighbors
          yn(k) = (spoints(1,1)*(k/radius)+origy+y(k))/2;
          xn(k) = (spoints(1,2)*(k/radius)+origx+x(k))/2;    
      else
          yn(k) = (spoints(i+1,1)*(k/radius)+origy+y(k))/2;
          xn(k) = (spoints(i+1,2)*(k/radius)+origx+x(k))/2;
      end
      if i==1
          yp(k) = (spoints(neighbors,1)*(k/radius)+origy+y(k))/2;
          xp(k) = (spoints(neighbors,2)*(k/radius)+origx+x(k))/2;  
      else
          yp(k) = (spoints(i-1,1)*(k/radius)+origy+y(k))/2;
          xp(k) = (spoints(i-1,2)*(k/radius)+origx+x(k))/2;
      end


      % Calculate floors, ceils and rounds for the x and y.
      fy(k) = floor(y(k)); cy(k) = ceil(y(k)); ry(k) = round(y(k));
      fx(k) = floor(x(k)); cx(k) = ceil(x(k)); rx(k) = round(x(k));
      fyn(k) = floor(yn(k)); cyn(k) = ceil(yn(k)); ry(k) = round(yn(k));
      fxn(k) = floor(xn(k)); cxn(k) = ceil(xn(k)); rx(k) = round(xn(k));
      fyp(k) = floor(yp(k)); cyp(k) = ceil(yp(k)); ryp(k) = round(yp(k));
      fxp(k) = floor(xp(k)); cxp(k) = ceil(xp(k)); rxp(k) = round(xp(k));
      % Check if interpolation is needed.
      if (abs(x(k) - rx(k)) < 1e-6) && (abs(y(k) - ry(k)) < 1e-6)
        % Interpolation is not needed, use original datatypes
        if k==1
            N(:,:,k) = image(ry(k):ry(k)+dy,rx(k):rx(k)+dx);
        else
            N(:,:,k) = (image(ry(k):ry(k)+dy,rx(k):rx(k)+dx)+image(ryn(k):ryn(k)+dy,rxn(k):rxn(k)+dx)+image(ryp(k):ryp(k)+dy,rxp(k):rxp(k)+dx))/3;
        end

      else
        % Interpolation needed, use double type images
        ty(k) = y(k) - fy(k);
        tx(k) = x(k) - fx(k);
        tyn(k) = yn(k) - fyn(k);
        txn(k) = xn(k) - fxn(k);
        typ(k) = yp(k) - fyp(k);
        txp(k) = xp(k) - fxp(k);

        % Compute interpolated pixel values
        a1(:,:,k)=d_image(fy(k):fy(k)+dy,fx(k):fx(k)+dx);
        a2(:,:,k)=d_image(fy(k):fy(k)+dy,cx(k):cx(k)+dx);
        a3(:,:,k)=d_image(cy(k):cy(k)+dy,fx(k):fx(k)+dx);
        a4(:,:,k)=d_image(cy(k):cy(k)+dy,cx(k):cx(k)+dx);

        a1n(:,:,k)=d_image(fyn(k):fyn(k)+dy,fxn(k):fxn(k)+dx);
        a2n(:,:,k)=d_image(fyn(k):fyn(k)+dy,cxn(k):cxn(k)+dx);
        a3n(:,:,k)=d_image(cyn(k):cyn(k)+dy,fxn(k):fxn(k)+dx);
        a4n(:,:,k)=d_image(cyn(k):cyn(k)+dy,cxn(k):cxn(k)+dx);

        a1p(:,:,k)=d_image(fyp(k):fyp(k)+dy,fxp(k):fxp(k)+dx);
        a2p(:,:,k)=d_image(fyp(k):fyp(k)+dy,cxp(k):cxp(k)+dx);
        a3p(:,:,k)=d_image(cyp(k):cyp(k)+dy,fxp(k):fxp(k)+dx);
        a4p(:,:,k)=d_image(cyp(k):cyp(k)+dy,cxp(k):cxp(k)+dx);

        horInterp1(:,:,k) =a1(:,:,k)*(1-tx(k))+a2(:,:,k)*tx(k);
        horInterp2(:,:,k) =a3(:,:,k)*(1-tx(k))+a4(:,:,k)*tx(k);
        horInterp1n(:,:,k) =a1n(:,:,k)*(1-txn(k))+a2n(:,:,k)*txn(k);
        horInterp2n(:,:,k) =a3n(:,:,k)*(1-txn(k))+a4n(:,:,k)*txn(k);
        horInterp1p(:,:,k) =a1p(:,:,k)*(1-txp(k))+a2p(:,:,k)*txp(k);
        horInterp2p(:,:,k) =a3p(:,:,k)*(1-txp(k))+a4p(:,:,k)*txp(k);
        if k==1
           N(:,:,k)=horInterp1(:,:,k)*(1-ty(k))+horInterp2(:,:,k)*ty(k);
        else
           N(:,:,k)=(horInterp1(:,:,k)*(1-ty(k))+horInterp2(:,:,k)*ty(k)+horInterp1n(:,:,k)*(1-tyn(k))+horInterp2n(:,:,k)*tyn(k)+horInterp1p(:,:,k)*(1-typ(k))+horInterp2p(:,:,k)*typ(k))/3; 
        end

      end
  end
    %Decorrelation
    if K==1
        Z1 = N(:,:,1);
        M1=abs(N(:,:,1)-d_C);
        A1=A1+Z1;
        MM1=MM1+M1;
    elseif K==2
        Z1=N(:,:,1)+N(:,:,2)-2*d_C ;
        Z2=N(:,:,1)-N(:,:,2);
        A1=A1+Z1;
        A2=A1+Z2;
        m1=abs(N(:,:,1)-d_C);
        m2=abs(N(:,:,2)-d_C);
        M1=m1+m2;
        M2=m1-m2;
        MM1=MM1+M1;
        MM2=MM2+M2;
    else
        Z1=N(:,:,1)+N(:,:,2)+N(:,:,3)-3*d_C;
        Z2=N(:,:,1)-N(:,:,3);
        Z3=N(:,:,1)-2*N(:,:,2)+N(:,:,3);
        A1=A1+Z1;
        A2=A2+Z2;
        A3=A2+Z3;
        m1=abs(N(:,:,1)-d_C);
        m2=abs(N(:,:,2)-d_C);
        m3=abs(N(:,:,3)-d_C);
        M1=m1+m2+m3;
        M2=m1-m2;
        M3=m1-2*m2+m3;
        MM1=MM1+M1;
        MM2=MM2+M2;
        MM3=MM3+M3;
    end

 end
%threshold
if K==1
    mA1=mean(mean(A1));
    mMM1=mean(mean(MM1));
elseif K==2
    mA1=mean(mean(A1));
    mA2=mean(mean(A2));
    mMM1=mean(mean(MM1));
    mMM2=mean(mean(MM2));
else
    mA1=mean(mean(A1));
    mA2=mean(mean(A2));
    mA3=mean(mean(A3));
    mMM1=mean(mean(MM1));
    mMM2=mean(mean(MM2));
    mMM3=mean(mean(MM3));
end
%Compute the dlbp
for i = 1:neighbors
  for k=1:radius  
  y(k) = spoints(i,1)*(k/radius)+origy;
  x(k) = spoints(i,2)*(k/radius)+origx;
  if i==neighbors
      yn(k) = (spoints(1,1)*(k/radius)+origy+y(k))/2;
      xn(k) = (spoints(1,2)*(k/radius)+origx+x(k))/2;    
  else
      yn(k) = (spoints(i+1,1)*(k/radius)+origy+y(k))/2;
      xn(k) = (spoints(i+1,2)*(k/radius)+origx+x(k))/2;
  end
  if i==1
      yp(k) = (spoints(neighbors,1)*(k/radius)+origy+y(k))/2;
      xp(k) = (spoints(neighbors,2)*(k/radius)+origx+x(k))/2;  
  else
      yp(k) = (spoints(i-1,1)*(k/radius)+origy+y(k))/2;
      xp(k) = (spoints(i-1,2)*(k/radius)+origx+x(k))/2;
  end
  
 
  % Calculate floors, ceils and rounds for the x and y.
  fy(k) = floor(y(k)); cy(k) = ceil(y(k)); ry(k) = round(y(k));
  fx(k) = floor(x(k)); cx(k) = ceil(x(k)); rx(k) = round(x(k));
  fyn(k) = floor(yn(k)); cyn(k) = ceil(yn(k)); ry(k) = round(yn(k));
  fxn(k) = floor(xn(k)); cxn(k) = ceil(xn(k)); rx(k) = round(xn(k));
  fyp(k) = floor(yp(k)); cyp(k) = ceil(yp(k)); ryp(k) = round(yp(k));
  fxp(k) = floor(xp(k)); cxp(k) = ceil(xp(k)); rxp(k) = round(xp(k));
  % Check if interpolation is needed.
  if (abs(x(k) - rx(k)) < 1e-6) && (abs(y(k) - ry(k)) < 1e-6)
    % Interpolation is not needed, use original datatypes
    if k==1
        N(:,:,k) = image(ry(k):ry(k)+dy,rx(k):rx(k)+dx);
    else
        N(:,:,k) = (image(ry(k):ry(k)+dy,rx(k):rx(k)+dx)+image(ryn(k):ryn(k)+dy,rxn(k):rxn(k)+dx)+image(ryp(k):ryp(k)+dy,rxp(k):rxp(k)+dx))/3;
    end
  else
    % Interpolation needed, use double type images
    ty(k) = y(k) - fy(k);
    tx(k) = x(k) - fx(k);
    tyn(k) = yn(k) - fyn(k);
    txn(k) = xn(k) - fxn(k);
    typ(k) = yp(k) - fyp(k);
    txp(k) = xp(k) - fxp(k);
    % Compute interpolated pixel values
    a1(:,:,k)=d_image(fy(k):fy(k)+dy,fx(k):fx(k)+dx);
    a2(:,:,k)=d_image(fy(k):fy(k)+dy,cx(k):cx(k)+dx);
    a3(:,:,k)=d_image(cy(k):cy(k)+dy,fx(k):fx(k)+dx);
    a4(:,:,k)=d_image(cy(k):cy(k)+dy,cx(k):cx(k)+dx);
    
    a1n(:,:,k)=d_image(fyn(k):fyn(k)+dy,fxn(k):fxn(k)+dx);
    a2n(:,:,k)=d_image(fyn(k):fyn(k)+dy,cxn(k):cxn(k)+dx);
    a3n(:,:,k)=d_image(cyn(k):cyn(k)+dy,fxn(k):fxn(k)+dx);
    a4n(:,:,k)=d_image(cyn(k):cyn(k)+dy,cxn(k):cxn(k)+dx);
    
    a1p(:,:,k)=d_image(fyp(k):fyp(k)+dy,fxp(k):fxp(k)+dx);
    a2p(:,:,k)=d_image(fyp(k):fyp(k)+dy,cxp(k):cxp(k)+dx);
    a3p(:,:,k)=d_image(cyp(k):cyp(k)+dy,fxp(k):fxp(k)+dx);
    a4p(:,:,k)=d_image(cyp(k):cyp(k)+dy,cxp(k):cxp(k)+dx);
    
    horInterp1(:,:,k) =a1(:,:,k)*(1-tx(k))+a2(:,:,k)*tx(k);
    horInterp2(:,:,k) =a3(:,:,k)*(1-tx(k))+a4(:,:,k)*tx(k);
    horInterp1n(:,:,k) =a1n(:,:,k)*(1-txn(k))+a2n(:,:,k)*txn(k);
    horInterp2n(:,:,k) =a3n(:,:,k)*(1-txn(k))+a4n(:,:,k)*txn(k);
    horInterp1p(:,:,k) =a1p(:,:,k)*(1-txp(k))+a2p(:,:,k)*txp(k);
    horInterp2p(:,:,k) =a3p(:,:,k)*(1-txp(k))+a4p(:,:,k)*txp(k);
    if k==1
       N(:,:,k)=horInterp1(:,:,k)*(1-ty(k))+horInterp2(:,:,k)*ty(k);
    else
       N(:,:,k)=(horInterp1(:,:,k)*(1-ty(k))+horInterp2(:,:,k)*ty(k)+horInterp1n(:,:,k)*(1-tyn(k))+horInterp2n(:,:,k)*tyn(k)+horInterp1p(:,:,k)*(1-typ(k))+horInterp2p(:,:,k)*typ(k))/3; 
    end
  end
  end
   v = 2^(i-1);
    if K==1
        Z1 = N(:,:,1) ;
        dLBP_D1=Z1>=mA1;
        resultd1 = resultd1 + v*dLBP_D1;
        M1=abs(N(:,:,1)-d_C);
        dLBP_M1=M1>=mMM1;
        resultM1 = resultM1 + v*dLBP_M1;
    elseif K==2
        Z1=N(:,:,1)+N(:,:,2)-2*d_C ;
        Z2=N(:,:,1)-N(:,:,2);
        dLBP_D1=Z1>=mA1;
        dLBP_D2=Z2>=mA2;
        resultd1 = resultd1 + v*dLBP_D1;
        resultd2 = resultd1 + v*dLBP_D2;
        m1=abs(N(:,:,1)-d_C);
        m2=abs(N(:,:,2)-d_C);
        M1=m1+m2;
        M2=m1-m2;
        dLBP_M1=M1>=mMM1;
        dLBP_M2=M2>=mMM2;
        resultM1 = resultM1 + v*dLBP_M1;
        resultM2 = resultM2 + v*dLBP_M2;
    else
        Z1=N(:,:,1)+N(:,:,2)+N(:,:,3)-3*d_C;
        Z2=N(:,:,1)-N(:,:,3);
        Z3=N(:,:,1)-2*N(:,:,2)+N(:,:,3);
        dLBP_D1=Z1>=mA1;
        dLBP_D2=Z2>=mA2;
        dLBP_D3=Z3>=mA3;
        resultd1 = resultd1 + v*dLBP_D1;
        resultd2 = resultd2 + v*dLBP_D2;
        resultd3 = resultd3 + v*dLBP_D3;
        m1=abs(N(:,:,1)-d_C);
        m2=abs(N(:,:,2)-d_C);
        m3=abs(N(:,:,3)-d_C);
        M1=m1+m2+m3;
        M2=m1-m2;
        M3=m1-2*m2+m3;
        dLBP_M1=M1>=mMM1;
        dLBP_M2=M2>=mMM2;
        dLBP_M3=M3>=mMM3;
        resultM1 = resultM1 + v*dLBP_M1;
        resultM2 = resultM2 + v*dLBP_M2;
        resultM3 = resultM3 + v*dLBP_M3;
    end

end
mmC=mean(mean(d_C));
resultC=d_C>=mmC;
%Apply mapping if it is defined
if isstruct(mapping)
    bins = mapping.num;
    for i = 1:size(resultd1,1)
        for j = 1:size(resultd1,2)
            if K==1
            resultd1(i,j) = mapping.table(resultd1(i,j)+1);
            resultM1(i,j) = mapping.table(resultM1(i,j)+1);
            elseif K==2
            resultd1(i,j) = mapping.table(resultd1(i,j)+1);
            resultM1(i,j) = mapping.table(resultM1(i,j)+1);
            resultd2(i,j) = mapping.table(resultd2(i,j)+1);
            resultM2(i,j) = mapping.table(resultM2(i,j)+1);
            else 
            resultd1(i,j) = mapping.table(resultd1(i,j)+1);
            resultM1(i,j) = mapping.table(resultM1(i,j)+1);
            resultd2(i,j) = mapping.table(resultd2(i,j)+1);
            resultM2(i,j) = mapping.table(resultM2(i,j)+1);
            resultd3(i,j) = mapping.table(resultd3(i,j)+1);
            resultM3(i,j) = mapping.table(resultM3(i,j)+1);
            end
        end
    end
end

if (strcmp(mode,'h') || strcmp(mode,'hist') || strcmp(mode,'nh'))
    % Return with LBP histogram if mode equals 'hist'.
    if K==1
         resultd1=hist(resultd1(:),0:(bins-1));
         resultM1=hist(resultM1(:),0:(bins-1));
         resultC=hist(resultC(:),0:1);
         result=[resultd1 resultM1 resultC];

    elseif K==2
         resultd1=hist(resultd1(:),0:(bins-1));
         resultM1=hist(resultM1(:),0:(bins-1));
         resultd2=hist(resultd2(:),0:(bins-1));
         resultM2=hist(resultM2(:),0:(bins-1));
         resultC=hist(resultC(:),0:1);
         result=[resultd1 resultM1 resultd2 resultM2 resultC];

    else
         resultd1=hist(resultd1(:),0:(bins-1));
         resultM1=hist(resultM1(:),0:(bins-1));
         resultd2=hist(resultd2(:),0:(bins-1));
         resultM2=hist(resultM2(:),0:(bins-1));
         resultd3=hist(resultd3(:),0:(bins-1));
         resultM3=hist(resultM3(:),0:(bins-1));
         resultC=hist(resultC(:),0:1);
         result=[resultd1 resultM1 resultd2 resultM2 resultd3 resultM3 resultC];

    end
    if (strcmp(mode,'nh'))
        result=result/sum(result);
    end
end

end


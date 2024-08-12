function [ result ] = ilbp( radius,n_sampling,img )
    [ysize, xsize] = size(img);
    I=img(1:ysize,1:xsize);
    A=func_edge_height(I);
    A=A(radius+1:ysize-radius,radius+1:xsize-radius);
    mapping=getmapping(n_sampling,'riu2');
    lbp=lbp2(I,radius,n_sampling,mapping,0);
    hlbp=zeros(1,mapping.num);
    hcss1=zeros(1,7);
    for i=0:mapping.num-1
        hlbp(1,i+1)=sum(A(lbp==i));
    end
    css1=css(I,radius,0);
    for i=0:6
        hcss1(1,i+1)=sum(A(css1==i));
    end
    result=[hlbp hcss1];
end

function result = css(varargin) % image,radius,neighbors,mapping,mode
narginchk(1,3);
radius=varargin{2};
img=varargin{1};
if(nargin >= 3)
    mode=varargin{3};
else
    mode='h';
end
s=radius;
Z = double(img);
[h, w, d] = size(Z);
result = zeros(h-2*s, w-2*s, d); % init intra-channel LBP images
T=2;
% obtain LBPs at every pixel for each intra-channel
for didx = 1:d
    

    X = zeros(6,h-2*s,w-2*s);
  I0= Z(1+s-s:h-s-s,1+s  :w-s,didx);
  I1=Z(1+s  :h-s  ,1+s-s:w-s-s,didx);
  I2=Z(1+s+s:h-s+s,1+s  :w-s,didx);
  I3=Z(1+s  :h-s  ,1+s+s:w-s+s,didx);

        X(1,:,:) = I2-I3;
        X(2,:,:) = I1-I3;
        X(3,:,:) = I1-I2;
        X(4,:,:) = I0-I3;
        X(5,:,:) = I0-I2;
        X(6,:,:) = I0-I1;

    X=double(X>T);
    result(:, :,didx) = reshape([1,2,4,8,16,32]*X(:,:), h-2*s, w-2*s);
 
end

%Apply mapping if it is defined
mapnum=getmapnum;
    bins = mapnum.num;
    for i = 1:size(result,1)
        for j = 1:size(result,2)
            result(i,j) = mapnum.table(result(i,j)+1);
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
%the CSPPs with same number of ¡®1¡¯ are more likely to present similar structure
function mapnum = getmapnum
n=2^6-1;
j=[0 0 0 0 0 0];
table=1:n+1;
for i=0:n
    table(i+1)=sum(bitget(bitxor(i,j),1:6));
end
mapnum.table=table;
mapnum.num=7;
end

%code from https://github.com/miaoxikui/JNDVP_for_IQA/blob/master/func_JND_modeling_pattern_complexity.m#L200
function edge_height = func_edge_height( img )
G1 = [0 0 0 0 0
   1 3 8 3 1
   0 0 0 0 0
   -1 -3 -8 -3 -1
   0 0 0 0 0];
G2=[0 0 1 0 0
   0 8 3 0 0
   1 3 0 -3 -1
   0 0 -3 -8 0
   0 0 -1 0 0];
G3=[0 0 1 0 0
   0 0 3 8  0
   -1 -3 0 3 1
   0 -8 -3 0 0
   0 0 -1 0 0];
G4=[0 1 0 -1 0
   0 3 0 -3 0
   0 8 0 -8 0
   0 3 0 -3 0
   0 1 0 -1 0];
% calculate the max grad
[size_x,size_y]=size(img);
grad=zeros(size_x,size_y,4);
grad(:,:,1) = filter2(G1,img)/16;
grad(:,:,2) = filter2(G2,img)/16;
grad(:,:,3) = filter2(G3,img)/16;
grad(:,:,4) = filter2(G4,img)/16;
max_gard = max( abs(grad), [], 3 );
maxgard = max_gard( 3:end-2, 3:end-2 );
edge_height =maxgard;
edge_height = padarray( maxgard, [2,2], 'symmetric' );
end
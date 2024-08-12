function [ result ] = lgbp(radius,n_sampling,img  )
img=imresize(img,[500,500]);
I=[img(:,:,1), img(:,:,2), img(:,:,3)];

g=load('data/gabor.mat','G');
%apply the filter to check the image
outMag = imgaborfilter(I,g);

%display the results
outSize = size(outMag);
mappingtype='riu2';
mode='h';
rowstep=32;
colstep=32;
for i=1:outSize(3)
    %lbpI=lbp2(outMag(:,:,i),radius,n_sampling,mapping,0);
    result1(i,:)= lbp_blk(outMag(:,:,i), radius, n_sampling, mappingtype,mode, rowstep, colstep);
end
result= result1(:);
result=result';
end


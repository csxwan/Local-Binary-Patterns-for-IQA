function [ result] = slbp( radius,n_sampling,img,smap  )
I=rgb2gray(img);
[x,y]=size(I);
mapping=0;
resultlbp=double(lbp2(I,radius,n_sampling,mapping,0));
bins = 2^n_sampling;
[xsize,ysize]=size(resultlbp);
% smap=[smap smap smap];
smap=double(smap(radius+1:x-radius,radius+1:y-radius));
result=zeros(1,bins);
for i=1:xsize
    for j=1:ysize
        result(1,resultlbp(i,j)+1)=result(1,resultlbp(i,j)+1)+smap(i,j);
    end
end
end


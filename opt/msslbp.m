function [ result] = msslbp( radius,n_sampling,img,smap  )
I=img;
[x,y]=size(I);
mapping=0;
radius1=1;
radius2=2;
radius3=3;
n_sampling4=4;
n_sampling8=8;
n_sampling16=16;
resultlbp14=double(lbp2(I,radius1,n_sampling4,mapping,0));
resultlbp18=double(lbp2(I,radius1,n_sampling8,mapping,0));
resultlbp28=double(lbp2(I,radius2,n_sampling8,mapping,0));
resultlbp216=double(lbp2(I,radius2,n_sampling16,mapping,0));
resultlbp38=double(lbp2(I,radius3,n_sampling8,mapping,0));
resultlbp316=double(lbp2(I,radius3,n_sampling16,mapping,0));
resultlbp24=double(lbp2(I,radius2,n_sampling4,mapping,0));
resultlbp34=double(lbp2(I,radius3,n_sampling4,mapping,0));
bins4 = 2^4;
bins8 = 2^8;
bins16 = 2^16;
[xsize1,ysize1]=size(resultlbp14);
smap=[smap smap smap];
smap1=double(smap(radius1+1:x-radius1,radius1+1:y-radius1));
smap2=double(smap(radius2+1:x-radius2,radius2+1:y-radius2));
smap3=double(smap(radius3+1:x-radius3,radius3+1:y-radius3));
result14=zeros(1,bins4);
result18=zeros(1,bins8);
result24=zeros(1,bins4);
result28=zeros(1,bins8);
result216=zeros(1,bins16);
result34=zeros(1,bins4);
result38=zeros(1,bins8);
result316=zeros(1,bins16);
for i=1:xsize1
    for j=1:ysize1
        result14(1,resultlbp14(i,j)+1)=result14(1,resultlbp14(i,j)+1)+smap1(i,j);
        result18(1,resultlbp18(i,j)+1)=result18(1,resultlbp18(i,j)+1)+smap1(i,j);
    end
end
[xsize2,ysize2]=size(resultlbp24);
for i=1:xsize2
    for j=1:ysize2
        result24(1,resultlbp24(i,j)+1)=result24(1,resultlbp24(i,j)+1)+smap2(i,j);
        result28(1,resultlbp28(i,j)+1)=result28(1,resultlbp28(i,j)+1)+smap2(i,j);
        result216(1,resultlbp216(i,j)+1)=result216(1,resultlbp216(i,j)+1)+smap2(i,j);
    end
end
[xsize3,ysize3]=size(resultlbp34);
for i=1:xsize3
    for j=1:ysize3
        result34(1,resultlbp34(i,j)+1)=result34(1,resultlbp34(i,j)+1)+smap3(i,j);
        result38(1,resultlbp38(i,j)+1)=result38(1,resultlbp38(i,j)+1)+smap3(i,j);
        result316(1,resultlbp316(i,j)+1)=result316(1,resultlbp316(i,j)+1)+smap3(i,j);
    end
end

result=[result14 result18 result24 result28 result216 result34 result38 result316];
end
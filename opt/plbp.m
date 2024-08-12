function [ result ] = plbp( radius,n_sampling,img )
I=img;
N=4;%pyramid levels
sigma=0.5;%Gaussian low pass filter with σ=0.5
windowsize=3;%Gaussian low pass filter with windowsize=3
subsample=2;
mapping=getmapping(n_sampling,'riu2');
result=[];
for i=1:N
    result=[result lbp2(I,radius,n_sampling,mapping,'h')];
    I=imresize(I,1/subsample);
    gausFilter = fspecial('gaussian',windowsize,sigma);      %matlab 自带高斯模板滤波
    I=imfilter(I,gausFilter,'replicate');
end


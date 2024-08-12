function [ resultp, resultn ] = ltp(radius,n_sampling,img  )
I=img;
mapping=getmapping(n_sampling,'riu2');
[resultp, resultn]=ltp2(I,radius,n_sampling,0,0,8);


end
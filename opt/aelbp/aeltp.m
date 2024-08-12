function [ result ] = aeltp( radius,n_sampling,img )
I=img;
mapping=getmapping(n_sampling,'riu2');
result = AE_ltp(I,radius,n_sampling,mapping.table,'h');  
end



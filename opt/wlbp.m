function [ result ] = wlbp( radius,n_sampling,img )
I=rgb2gray(img);
[LL,LH,HL,HH] = dwt2(I,'haar');
 mapping=getmapping(n_sampling,'ri');
 %WLBP-I
 resultHL = lbpT( HL,radius, n_sampling, mapping,'nh',1);
 resultHH = lbpT( HH,radius, n_sampling, mapping,'nh',1);
%  %WLBP-II
%  resultLH = lbpT( LH,radius, n_sampling, mapping,'nh',1);
%  resultHL = lbpT( HL,radius, n_sampling, mapping,'nh',1);
%  resultHH = lbpT( HLH,radius, n_sampling, mapping,'nh',1);
%  %WLBP-III
  resultLL = lbpT( LL,radius, n_sampling, mapping,'nh',4);
%  resultLH = lbpT( LH,radius, n_sampling, mapping,'nh',1);
%  resultHL = lbpT( HL,radius, n_sampling, mapping,'nh',1);
%  resultHH = lbpT( HH,radius, n_sampling, mapping,'nh',1);
result=[resultLL resultHL resultHH];
end


function [ result ] = mslbp(~,~,img ,mapping4 ,mapping8,mapping16,mapping24)
I=img;

result1=lbp2(I,1,4,mapping4,'h');
result2=lbp2(I,1,8,mapping8,'h');
result3=lbp2(I,2,4,mapping4,'h');
result4=lbp2(I,2,8,mapping8,'h');
result5=lbp2(I,2,16,mapping16,'h');
result6=lbp2(I,3,4,mapping4,'h');
result7=lbp2(I,3,8,mapping8,'h');
result8=lbp2(I,3,16,mapping16,'h');
result9=lbp2(I,3,24,mapping24,'h');
result=[result1 result2 result3 result4 result5 result6 result7 result8 result9];
end


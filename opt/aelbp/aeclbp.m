function [ result ] = aeclbp( radius,n_sampling,img )
I=img;
mapping=getmapping(n_sampling,'riu2');
[AECLBP_S,AECLBP_M,AECLBP_C] = AE_clbp(I,radius,n_sampling,mapping,'x');
    
% Generate histogram of AECLBP_S
AECLBP_SH = hist(AECLBP_S(:),0:mapping.num-1);

% Generate histogram of AECLBP_M
AECLBP_MH = hist(AECLBP_M(:),0:mapping.num-1);    

% Generate histogram of AECLBP_M/C
AECLBP_MC = [AECLBP_M(:),AECLBP_C(:)];
Hist3D = hist3(AECLBP_MC,[mapping.num,2]);
AECLBP_MCH = reshape(Hist3D,1,numel(Hist3D));

% Generate histogram of AECLBP_S_M/C
AECLBP_S_MCH = [AECLBP_SH,AECLBP_MCH];

% Generate histogram of AECLBP_S/M
AECLBP_SM = [AECLBP_S(:),AECLBP_M(:)];
Hist3D = hist3(AECLBP_SM,[mapping.num,mapping.num]);
AECLBP_SMH = reshape(Hist3D,1,numel(Hist3D));

% Generate histogram of AECLBP_S/M/C
AECLBP_MCSum = AECLBP_M;
idx = find(AECLBP_C);
AECLBP_MCSum(idx) = AECLBP_MCSum(idx)+mapping.num;
AECLBP_SMC = [AECLBP_S(:),AECLBP_MCSum(:)];
Hist3D = hist3(AECLBP_SMC,[mapping.num,mapping.num*2]);
AECLBP_SMCH = reshape(Hist3D,1,numel(Hist3D));
result=AECLBP_SMCH;

end



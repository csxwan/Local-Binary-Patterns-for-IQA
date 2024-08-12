 function dataset = getDataset( datasetName )
% datasetName can be one of the followings:
%     'TID2008'       
%     'TID2013'
%     'CSIQ'
%     'LIVE'
%     'IVC'
%     'MICT'
%     'WIQ'
%     'A57'
%     'LIVE_C'
%
% fields name of dataset are
%   img_path_ref
%   img_path_dist
%   to_ref
%   scores
%   scores_type

dataset = feval( [ 'get' datasetName ] );

end


function dataset = getA57()
lib_path='dataset/a57';
load dataset/a57.mat

lib_ref=strcat(lib_path,'/','src_imgs');
src={'horse', 'harbour','baby'};
db_deg={'A', 'B','C'};
imgs_type={'FLT', 'JPG','JP2','DCQ','BLR','NOZ'};
num_ref=length(src);
lib_dist=strcat(lib_path,'/','dst_imgs'); % library of distortion images

num_dist = num_ref*length(db_deg)*length(imgs_type);

images_path_ref = cell(num_ref,1);
images_path_dist = cell(num_dist,1);
to_ref = zeros(num_dist,1);

cnt_dist = 0;
for idx_ref = 1:num_ref
    images_path_ref{idx_ref}=strcat(lib_ref,'/',src{idx_ref},'.bmp');
    for idx_deg=1:length(db_deg)
        for idx_type=1:length(imgs_type)
            cnt_dist = cnt_dist + 1;
            to_ref(cnt_dist) = idx_ref;
            images_path_dist{cnt_dist}=strcat(lib_dist,'/',src{idx_ref},'/',db_deg{idx_deg},'/',src{idx_ref},'.',imgs_type{idx_type},'.bmp');
        end
    end
end

dataset.img_path_ref = images_path_ref;
dataset.img_path_dist = images_path_dist;
dataset.to_ref = to_ref;
dataset.scores = a57.mos;
dataset.scores_type = 'mos'; % ???

end % function getA57


function dataset = getWIQ()
lib_path='dataset/wiq';
load dataset/wiq_info.mat

dist_num =[5 6 9 24 11 12 13];
% idx_ref=[4,4,4,4,4,5,4,4,5,4,4,7,3,6,5,2,4,2,6,4,...
%     4,7,4,2,4,5,4,6,7,3,6,6,4,4,2,4,4,7,3,4,...
%     3,7,7,6,5,7,7,1,5,6,1,5,2,2,7,6,3,7,5,6,...
%     5,7,6,7,6,7,1,3,1,3,4,6,4,5,5,1,4,3,3,4];

num_ref = 7;
num_dist = 80;

images_path_ref = cell(num_ref,1);
images_path_dist = cell(num_dist,1);
to_ref = zeros(num_dist,1);

cnt_dist = 0;
for idx_ref=1:num_ref
    images_path_ref{idx_ref}=strcat( lib_path,'/','ref_img_',num2str(idx_ref,'%.3d'),'.bmp');
    cum_sum=sum(dist_num(1:idx_ref));
    sub_start=cum_sum-dist_num(idx_ref)+1;
    sub_end=cum_sum;
    for order_dist=sub_start:sub_end;
        cnt_dist = cnt_dist + 1;
        to_ref(cnt_dist) = idx_ref;
        images_path_dist{cnt_dist}=strcat(lib_path,wiq_info.img_dist{order_dist}) ;
    end

end

dataset.img_path_ref = images_path_ref;
dataset.img_path_dist = images_path_dist;
dataset.to_ref = to_ref;
dataset.scores = wiq_info.mos;
dataset.scores_type = 'mos';

end % function getWIQ


function dataset = getMICT()
lib_path='dataset/mict';
mict=load('dataset/mict_dist.mat');
mict = mict.mict;

num_ref = 28;
num_dist = length(mict); % 168

images_path_ref = cell(num_ref,1);
images_path_dist = cell(num_dist,1);
to_ref = zeros(num_dist,1);

cnt_ref=0;
tmp_pre=[];
for idx_dist=1:num_dist
    tmp=mict(idx_dist).reference_image_name ;
    if ~strcmp(tmp_pre,tmp )
        cnt_ref=cnt_ref+1;
        if idx_dist<=(num_dist/2) % 84
            subdir = 'jpeg';
        else
            subdir = 'jpeg2000';
        end
        images_path_ref{cnt_ref}=[ lib_path, filesep, subdir, filesep, tmp ];
        tmp_pre=tmp;
    end
end % loop for find images_path_ref

cnt_dist = 0;
for idx_ref=1:num_ref
    cum_sum=6*idx_ref;
    sub_start=cum_sum-6+1;
    sub_end=cum_sum;
    for order_dist=sub_start:sub_end;
        cnt_dist = cnt_dist + 1;
        to_ref(cnt_dist) = idx_ref;
        if order_dist<=(num_dist/2) % 84
            subdir = 'jpeg';
        else
            subdir = 'jpeg2000';
        end
        images_path_dist{cnt_dist}=[ lib_path, filesep, subdir, filesep, mict(order_dist).distored_image_name ];
    end
end % loop for find images_path_dist

dataset.img_path_ref = images_path_ref;
dataset.img_path_dist = images_path_dist;
dataset.to_ref = to_ref;
dataset.scores = cell2mat( {mict.mos} )';
dataset.scores_type = 'mos';

end % function getMICT


function dataset = getIVC()
lib_path='dataset/ivc/color';
load 'dataset/ivc'

lib_ref=lib_path;%strcat(lib_path,'/','color');
src={'avion', 'barba','boats','clown','fruit','house','isabe','lenat','mandr','pimen'};
imgs_type={'j2000','jpeg','flou','jpeg_lumichr','lar'};

num_ref=length(src);
num_dist = 185;

images_path_ref = cell(num_ref,1);
images_path_dist = cell(num_dist,1);
to_ref = zeros(num_dist,1);

cnt_dist = 0;
for idx_ref=1:num_ref;
    images_path_ref{idx_ref}=strcat(lib_ref,'/',src{idx_ref},'.bmp');
    for idx_type=1:length(imgs_type);
        for idx_deg=1:5
            if strcmp(imgs_type{idx_type},'flou')
                image_name=strcat(lib_ref,'/',src{idx_ref},'_',imgs_type{idx_type},'_f',num2str(idx_deg),'.bmp');
            else
                image_name=strcat(lib_ref,'/',src{idx_ref},'_',imgs_type{idx_type},'_r',num2str(idx_deg),'.bmp');
            end
            if ~exist(image_name,'file')
                continue;
            end
            cnt_dist = cnt_dist + 1;
            to_ref(cnt_dist) = idx_ref;
            images_path_dist{cnt_dist}=image_name;
        end %for distored degree
    end %for distored
end %for original

dataset.img_path_ref = images_path_ref;
dataset.img_path_dist = images_path_dist;
dataset.to_ref = to_ref;
dataset.scores = ivc.mos;
dataset.scores_type = 'mos';

end % function getIVC


function dataset = getLIVE()
lib_path='dataset/live';
% live = load( 'dataset/live_data.mat' );
% live = live.live;
dmos = load( 'dataset/live/dmos.mat');
% orgs=dmos.orgs;
dmos=dmos.dmos;
refnames_all=load('dataset/live/refnames_all.mat');
refnames_all=refnames_all.refnames_all;
num_ref = 29;
num_dist = 982;

images_path_ref = cell(num_ref,1);
images_path_dist = cell(num_dist,1);
to_ref = zeros(num_dist,1);
dmos1=zeros(num_dist,1);
% imh = cell(num_dist,1);
% lib=dir(lib_path);%information of lib directory
% tlib_num=length(lib); %num of lib files

% cnt_ref = 0;
cnt_dist = 0;

 dist_name={'jp2k','jpeg','wn','gblur','fastfading'};
 dist_num=[227,233,174,174,174];
 for cnt_ref = 1:num_ref;
    images_ref=dir(fullfile('dataset/live/refimgs','/*.bmp'));
%   len=length(images); %images in this lib
    images_path_ref{cnt_ref}=strcat('dataset/live/refimgs','/',images_ref(cnt_ref).name); 
    name=images_ref(cnt_ref).name;
 for lib_num=1:5 % exclude '.' and '..'
     lib_name=strcat(lib_path,'/',dist_name{lib_num});
%     if lib_name==0  %whether the lib is exist
%         disp('file not exit!');
%     end
    images=dir(fullfile(lib_name,'/*.bmp'));
    DirCell=struct2cell(images);%按照数字进行排序
    imagess = sort_nat(DirCell(1,:));
%     imagess()=sort_nat({images().name});
%      len=length(images); %images in this lib
%      for i=1:len
%      imh{i}=imagess(i);
%      end
    len=dist_num(lib_num);
    cum_sum=sum(dist_num(1:lib_num));
    sub_start=cum_sum- dist_num(lib_num)+1;
%          sub_end=cum_sum;

    for  img_num=1:len 
        order=sub_start+img_num-1;
        if  strcmp(refnames_all{order},name)
%        if  strcmp(ref_name{order},name)
           cnt_dist = cnt_dist + 1;
           to_ref(cnt_dist) = cnt_ref;
           images_path_dist{cnt_dist}=strcat(lib_name,'/',imagess{img_num}); 
           dmos1(cnt_dist) = dmos(order);
%             imh{cnt_dist}=imagess(img_num).name;
           
       end
    end    
  
 end
end
dataset.img_path_ref = images_path_ref;
dataset.img_path_dist = images_path_dist;
   save('dataset\images_path_ref.mat','images_path_ref');
    save('dataset\images_path_dist.mat','images_path_dist');
dataset.to_ref = to_ref;
    save('dataset\to_ref.mat','to_ref');
%    save('dataset\imh.mat','imh');

dataset.scores = dmos1;
  save('dataset\dmos1.mat','dmos1');
dataset.scores_type = 'dmos'; % ???

end % function getLIVE
%sort_nat具体内容

function [cs,index] = sort_nat(c,mode)
%sort_nat: Natural order sort of cell array of strings.
% usage:  [S,INDEX] = sort_nat(C)
%
% where,
%    C is a cell array (vector) of strings to be sorted.
%    S is C, sorted in natural order.
%    INDEX is the sort order such that S = C(INDEX);
%
% Natural order sorting sorts strings containing digits in a way such that
% the numerical value of the digits is taken into account.  It is
% especially useful for sorting file names containing index numbers with
% different numbers of digits.  Often, people will use leading zeros to get
% the right sort order, but with this function you don't have to do that.
% For example, if C = {'file1.txt','file2.txt','file10.txt'}, a normal sort
% will give you
%
%       {'file1.txt'  'file10.txt'  'file2.txt'}
%
% whereas, sort_nat will give you
%
%       {'file1.txt'  'file2.txt'  'file10.txt'}
%
% See also: sort

% Version: 1.4, 22 January 2011
% Author:  Douglas M. Schwarz
% Email:   dmschwarz=ieee*org, dmschwarz=urgrad*rochester*edu
% Real_email = regexprep(Email,{'=','*'},{'@','.'})


% Set default value for mode if necessary.
if nargin < 2
    mode = 'ascend';
end

% Make sure mode is either 'ascend' or 'descend'.
modes = strcmpi(mode,{'ascend','descend'});
is_descend = modes(2);
if ~any(modes)
    error('sort_nat:sortDirection',...
        'sorting direction must be ''ascend'' or ''descend''.')
end

% Replace runs of digits with '0'.
c2 = regexprep(c,'\d+','0');

% Compute char version of c2 and locations of zeros.
s1 = char(c2);
z = s1 == '0';

% Extract the runs of digits and their start and end indices.
[digruns,first,last] = regexp(c,'\d+','match','start','end');

% Create matrix of numerical values of runs of digits and a matrix of the
% number of digits in each run.
num_str = length(c);
max_len = size(s1,2);
num_val = NaN(num_str,max_len);
num_dig = NaN(num_str,max_len);
for i = 1:num_str
    num_val(i,z(i,:)) = sscanf(sprintf('%s ',digruns{i}{:}),'%f');
    num_dig(i,z(i,:)) = last{i} - first{i} + 1;
end

% Find columns that have at least one non-NaN.  Make sure activecols is a
% 1-by-n vector even if n = 0.
activecols = reshape(find(~all(isnan(num_val))),1,[]);
n = length(activecols);

% Compute which columns in the composite matrix get the numbers.
numcols = activecols + (1:2:2*n);

% Compute which columns in the composite matrix get the number of digits.
ndigcols = numcols + 1;

% Compute which columns in the composite matrix get chars.
charcols = true(1,max_len + 2*n);
charcols(numcols) = false;
charcols(ndigcols) = false;

% Create and fill composite matrix, comp.
comp = zeros(num_str,max_len + 2*n);
comp(:,charcols) = double(s1);
comp(:,numcols) = num_val(:,activecols);
comp(:,ndigcols) = num_dig(:,activecols);

% Sort rows of composite matrix and use index to sort c in ascending or
% descending order, depending on mode.
[~,index] = sortrows(comp);
if is_descend
    index = index(end:-1:1);
end
index = reshape(index,size(c));
cs = c(index);
end


function dataset = getCSIQ()
lib_path='dataset/csiq';
csiq=load ('dataset/csiq_dist.mat');
dmos=load( 'dataset/csiq.mat' );
dmos=dmos.csiq.dmos;

img_ref={ '1600.png' 'aerial_city.png' 'boston.png' 'bridge.png' 'butter_flower.png' 'cactus.png' 'child_swimming.png' 'couple.png' 'elk.png'...
    'family.png' 'fisher.png' 'foxy.png' 'geckos.png' 'lady_liberty.png' 'lake.png' 'log_seaside.png' 'monument.png' 'native_american.png'...
    'redwood.png' 'roping.png' 'rushmore.png' 'shroom.png' 'snow_leaves.png' 'sunset_sparrow.png' 'sunsetcolor.png' 'swarm.png' 'trolley.png'...
    'turtle.png' 'veggies.png' 'woman.png'};
dist_num=[ 29 29 29 29 29 29 29 29 28 29 29 29 28 28 29 29 29 29 29 29 29 29 29 29 29 28 29 29 29 29];

d = dir( [ lib_path, filesep, 'distored' ] );
dist_sub_dirs = { d(3:end).name };
dist_sub_dirs = sort( dist_sub_dirs ); % ascend, since sort does not support MODE
dist_sub_dirs = dist_sub_dirs(end:-1:1); % descend

num_ref=length(img_ref);
num_dist = 866;

images_path_ref = cell(num_ref,1);
images_path_dist = cell(num_dist,1);
to_ref = zeros(num_dist,1);

img_dist=csiq.img_dist;

cnt_dist = 0;
for idx_ref=1:num_ref
    images_path_ref{idx_ref}=[ lib_path, filesep, 'reference', filesep, img_ref{idx_ref} ];
    
    cum_sum=sum(dist_num(1:idx_ref));
    sub_start=cum_sum- dist_num(idx_ref)+1;
    sub_end=cum_sum;
    for order_dist=sub_start:sub_end;
        name1 = img_dist(order_dist).name_uc;
        name2 = img_dist(order_dist).name;
        
        % find sub directory
        sub_dir = '';
        for i = 1:numel(dist_sub_dirs)
            % jpeg2000 should be test before jpeg, so dist_sub_dirs should be descend
            if ~isempty( strfind( lower([name1,name2]),lower(dist_sub_dirs{i}) ) )
                sub_dir = dist_sub_dirs{i};
                break;
            end
        end
        
        % find correct path
        dir_path = [ lib_path, filesep, 'distored', filesep, sub_dir ];
        img_path_dist = [ dir_path, filesep, name1 ];
        if ~exist(img_path_dist,'file')
            img_path_dist = [ dir_path, filesep, name2 ];
        end

        cnt_dist = cnt_dist + 1;
        to_ref(cnt_dist) = idx_ref;
        images_path_dist{cnt_dist}=img_path_dist;
    end
end % for idx_ref

dataset.img_path_ref = images_path_ref;
dataset.img_path_dist = images_path_dist;
dataset.to_ref = to_ref;
dataset.scores = dmos;
dataset.scores_type = 'dmos';
  save('dataset\csiq\images_path_ref.mat','images_path_ref');
  save('dataset\csiq\images_path_dist.mat','images_path_dist');
  save('dataset\csiq\to_ref.mat','to_ref');
  save('dataset\csiq\dmos.mat','dmos');
end % function getCSIQ


function dataset = getTID2008()
lib_path='dataset/tid2008';
load dataset/tid2008/tid2008.mat;

lib_ref=strcat(lib_path,'/','reference_images');
lib_dist=strcat(lib_path,'/','distorted_images');

imgs_ref=dir(fullfile(lib_ref));

num_ref=length(imgs_ref)-2; % 25
num_dist = 1700;

images_path_ref = cell(num_ref,1);
images_path_dist = cell(num_dist,1);
to_ref = zeros(num_dist,1);

cnt_dist = 0;
for idx_ref=1:num_ref
    image_path=strcat(lib_ref,'/','I',num2str(idx_ref,'%.2d'),'.BMP');
    if ~exist(image_path,'file')
        image_path=strcat(lib_ref,'/','i',num2str(idx_ref,'%.2d'),'.bmp');
    end
    images_path_ref{idx_ref}=image_path;
    
    for idx_type=1:17
        for idx_deg=1:4
            image_name=strcat(lib_dist,'/','I',num2str(idx_ref,'%.2d'),'_',num2str(idx_type,'%.2d'),'_',num2str(idx_deg),'.bmp');
            if ~exist(image_name,'file')
                image_name=strcat(lib_dist,'/','i',num2str(idx_ref,'%.2d'),'_',num2str(idx_type,'%.2d'),'_',num2str(idx_deg),'.bmp');
            end
            
            if ~exist(image_name,'file')
                disp('error reading image file !');
            end
            
            cnt_dist = cnt_dist + 1;
            to_ref(cnt_dist) = idx_ref;
            images_path_dist{cnt_dist}=image_name;
        end %for distored degree
    end %for distored

end %for original

dataset.img_path_ref = images_path_ref;
dataset.img_path_dist = images_path_dist;
dataset.to_ref = to_ref;
dataset.scores = mos;
dataset.scores_type = 'mos';

end % function getTID2008


function dataset = getTID2013()
lib_path='dataset/tid2013';
load dataset/tid2013/tid2013.mat;
%load dataset/tid2013/mos.mat
%load ('E:\Users\a\temp\SRP\IQA_chi_nodatasets\dataset\tid2013\mos.mat')


num_ref=25;
num_dist = num_ref*24*5; % 3000

images_path_ref = cell(num_ref,1);
images_path_dist = cell(num_dist,1);
to_ref = zeros(num_dist,1);

cnt_dist = 0;
for idx_ref=1:num_ref
    image_path=strcat(lib_path,'/reference_images','/','I',num2str(idx_ref,'%.2d'),'.BMP');
    if ~exist(image_path,'file')
        image_path=strcat(lib_path,'/reference_images','/','i',num2str(idx_ref,'%.2d'),'.bmp');
    end
    images_path_ref{idx_ref} = image_path;

    for idx_type=1:24
        for idx_deg=1:5
            image_name=strcat(lib_path,'/distorted_images','/','I',num2str(idx_ref,'%.2d'),'_',num2str(idx_type,'%.2d'),'_',num2str(idx_deg),'.bmp');
            if ~exist(image_name,'file')
                image_name=strcat(lib_path,'/distorted_images','/','i',num2str(idx_ref,'%.2d'),'_',num2str(idx_type,'%.2d'),'_',num2str(idx_deg),'.bmp');
            end
            if ~exist(image_name,'file')
                image_name=strcat(lib_path,'/distorted_images','/','I',num2str(idx_ref,'%.2d'),'_',num2str(idx_type,'%.2d'),'_',num2str(idx_deg),'.BMP');
            end
            if ~exist(image_name,'file')
                image_name=strcat(lib_path,'/distorted_images','/','i',num2str(idx_ref,'%.2d'),'_',num2str(idx_type,'%.2d'),'_',num2str(idx_deg),'.BMP');
            end
            
            if ~exist(image_name,'file')
                disp('error reading image file !');
            end
            cnt_dist = cnt_dist + 1;
            to_ref(cnt_dist) = idx_ref;
            images_path_dist{cnt_dist}=image_name; 
        end
    end
end

dataset.img_path_ref = images_path_ref;
dataset.img_path_dist = images_path_dist;
dataset.to_ref = to_ref;
save('dataset\tid2013\to_ref.mat','to_ref');
dataset.scores = mos;
dataset.scores_type = 'mos';

end % function getTID2013


function dataset = getLIVE_C()
% lib_path='dataset/live_c';
mos = load( 'dataset/live_c/mos.mat');
mos=mos.mos2;
images_dist=dir(fullfile('dataset/live_c/Images'));
num_dist = length(images_dist);
images_path_dist = cell(num_dist,1);
% to_ref = zeros(num_dist,1);
% cnt_dist = 0;
% DirCell=struct2cell(images_dist);%按照数字进行排序
% images_dist = sort_nat(DirCell(1,:));
 for cnt_dist = 3:num_dist;
%     images_dist=dir(fullfile('dataset/live_c/Images','/*.bmp'));
%     images_dist=dir(fullfile('dataset/live_c/Images'));
%     images_dist=images_dist(3:end);
%   len=length(images); %images in this lib
    images_path_dist{cnt_dist}=strcat('dataset/live_c/Images','/',images_dist(cnt_dist).name); 

 end
 images_path_dist=images_path_dist(3:end);
dataset.img_path_dist = images_path_dist;
    save('dataset\images_path_dist.mat','images_path_dist');
dataset.scores = mos;
%   save('dataset\dmos1.mat','dmos1');
dataset.scores_type = 'mos'; % ???
dataset.img_path_ref = [];
dataset.to_ref = [];
end % function getLIVE_C


function dataset = getKonIQ_10k()
% lib_path='dataset/live_c';
mos = load( 'dataset/KonIQ-10k/mos.mat');
mos=mos.MOS;
images_dist=dir(fullfile('dataset/KonIQ-10k/50'));
num_dist = length(images_dist);
images_path_dist = cell(num_dist,1);
% to_ref = zeros(num_dist,1);
% cnt_dist = 0;
% DirCell=struct2cell(images_dist);%°?????×????????ò
% images_dist = sort_nat(DirCell(1,:));
 for cnt_dist = 3:num_dist
%     images_dist=dir(fullfile('dataset/live_c/Images','/*.bmp'));
%     images_dist=dir(fullfile('dataset/live_c/Images'));
%     images_dist=images_dist(3:end);
%   len=length(images); %images in this lib
    images_path_dist{cnt_dist}=strcat('dataset/KonIQ-10k/50','/',images_dist(cnt_dist).name); 

 end
 images_path_dist=images_path_dist(3:end);
dataset.img_path_dist = images_path_dist;
    save('dataset\images_path_dist.mat','images_path_dist');
% dataset.scores = mos;
dataset.scores = [];
%   save('dataset\dmos1.mat','dmos1');
dataset.scores_type = 'mos'; % ???
dataset.img_path_ref = [];
dataset.to_ref = [];
end % function getKonIQ-10k

function dataset = getSPAQ()
mos = load( 'dataset/SPAQ/SPAQmos.mat');
mos=mos.MOS;
images=dir(fullfile('dataset/SPAQ/testimages'));
DirCell=struct2cell(images);%按照数字进行排序
images_dist = sort_nat(DirCell(1,:));
num_dist = length(images_dist);
 for cnt_dist = 3:num_dist
    images_path_dist{cnt_dist}=strcat('dataset/SPAQ/testimages','/',images_dist(cnt_dist).name); 
 end
 images_path_dist=images_path_dist(3:end);
dataset.img_path_dist = images_path_dist;
save('dataset\images_path_dist.mat','images_path_dist');
dataset.scores =mos;
dataset.scores_type = 'mos'; 
dataset.img_path_ref = [];
dataset.to_ref = [];

end % function getSPAQ


function dataset = getblur()
lib_path='dataset/csiq';
csiq=load ('dataset/csiq_dist.mat');
dmos=load( 'dataset/csiq.mat' );
label=dmos.csiq.label;
dmos=dmos.csiq.dmos;


img_ref={ '1600.png' 'aerial_city.png' 'boston.png' 'bridge.png' 'butter_flower.png' 'cactus.png' 'child_swimming.png' 'couple.png' 'elk.png'...
    'family.png' 'fisher.png' 'foxy.png' 'geckos.png' 'lady_liberty.png' 'lake.png' 'log_seaside.png' 'monument.png' 'native_american.png'...
    'redwood.png' 'roping.png' 'rushmore.png' 'shroom.png' 'snow_leaves.png' 'sunset_sparrow.png' 'sunsetcolor.png' 'swarm.png' 'trolley.png'...
    'turtle.png' 'veggies.png' 'woman.png'};
% dist_num=[5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5];

num_ref=length(img_ref);
num_dist = 150;

images_path_ref = cell(num_ref,1);
images_path_dist = cell(num_dist,1);
to_ref = zeros(num_dist,1);
mos = zeros(num_dist,1);

img_dist=csiq.img_dist;


cnt_dist = 0;
for idx_ref=1:num_ref
    images_path_ref{idx_ref}=[ lib_path, filesep, 'reference', filesep, img_ref{idx_ref} ];
    
%     cum_sum=sum(dist_num(1:idx_ref));
%     sub_start=cum_sum- dist_num(idx_ref)+1;
%     sub_end=cum_sum;
%     for order_dist=sub_start:sub_end;
     for order_dist=1:866;
        name1 = img_dist(order_dist).name_uc;
        str=img_ref{idx_ref};
        len=size(str);
        str1=name1(1:len(2)-5);
        str2=str(1:len(2)-5);
        if(label(order_dist)==5 &&strcmp(str1,str2)==1)
%         name1 = img_dist(order_dist).name_uc;
        name2 = img_dist(order_dist).name;
        
        % find sub directory
        sub_dir = 'blur';
        
        % find correct path
        dir_path = [ lib_path, filesep, 'distored', filesep, sub_dir ];
        img_path_dist = [ dir_path, filesep, name1 ];
        if ~exist(img_path_dist,'file')
            img_path_dist = [ dir_path, filesep, name2 ];
        end

        cnt_dist = cnt_dist + 1;
        to_ref(cnt_dist) = idx_ref;
        images_path_dist{cnt_dist}=img_path_dist;
        mos(cnt_dist)=dmos(order_dist);
        end
    end
end % for idx_ref

dataset.img_path_ref = images_path_ref;
dataset.img_path_dist = images_path_dist;
dataset.to_ref = to_ref;
dataset.scores = mos;
dataset.scores_type = 'dmos';
end % function getCSIQ

function dataset = getawgn()
lib_path='dataset/csiq';
csiq=load ('dataset/csiq_dist.mat');
dmos=load( 'dataset/csiq.mat' );
label=dmos.csiq.label;
dmos=dmos.csiq.dmos;


img_ref={ '1600.png' 'aerial_city.png' 'boston.png' 'bridge.png' 'butter_flower.png' 'cactus.png' 'child_swimming.png' 'couple.png' 'elk.png'...
    'family.png' 'fisher.png' 'foxy.png' 'geckos.png' 'lady_liberty.png' 'lake.png' 'log_seaside.png' 'monument.png' 'native_american.png'...
    'redwood.png' 'roping.png' 'rushmore.png' 'shroom.png' 'snow_leaves.png' 'sunset_sparrow.png' 'sunsetcolor.png' 'swarm.png' 'trolley.png'...
    'turtle.png' 'veggies.png' 'woman.png'};
% dist_num=[5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5];

num_ref=length(img_ref);
num_dist = 150;

images_path_ref = cell(num_ref,1);
images_path_dist = cell(num_dist,1);
to_ref = zeros(num_dist,1);
mos = zeros(num_dist,1);

img_dist=csiq.img_dist;


cnt_dist = 0;
for idx_ref=1:num_ref
    images_path_ref{idx_ref}=[ lib_path, filesep, 'reference', filesep, img_ref{idx_ref} ];
    
%     cum_sum=sum(dist_num(1:idx_ref));
%     sub_start=cum_sum- dist_num(idx_ref)+1;
%     sub_end=cum_sum;
%     for order_dist=sub_start:sub_end;
     for order_dist=1:866;
        name1 = img_dist(order_dist).name_uc;
        str=img_ref{idx_ref};
        len=size(str);
        str1=name1(1:len(2)-5);
        str2=str(1:len(2)-5);
        if(label(order_dist)==1 &&strcmp(str1,str2)==1)
%         name1 = img_dist(order_dist).name_uc;
        name2 = img_dist(order_dist).name;
        
        % find sub directory
        sub_dir = 'awgn';
        
        % find correct path
        dir_path = [ lib_path, filesep, 'distored', filesep, sub_dir ];
        img_path_dist = [ dir_path, filesep, name1 ];
        if ~exist(img_path_dist,'file')
            img_path_dist = [ dir_path, filesep, name2 ];
        end

        cnt_dist = cnt_dist + 1;
        to_ref(cnt_dist) = idx_ref;
        images_path_dist{cnt_dist}=img_path_dist;
        mos(cnt_dist)=dmos(order_dist);
        end
    end
end % for idx_ref

dataset.img_path_ref = images_path_ref;
dataset.img_path_dist = images_path_dist;
dataset.to_ref = to_ref;
dataset.scores = mos;
dataset.scores_type = 'dmos';
end % function getCSIQ


function dataset = getcontrast()
lib_path='dataset/csiq';
csiq=load ('dataset/csiq_dist.mat');
dmos=load( 'dataset/csiq.mat' );
label=dmos.csiq.label;
dmos=dmos.csiq.dmos;


img_ref={ '1600.png' 'aerial_city.png' 'boston.png' 'bridge.png' 'butter_flower.png' 'cactus.png' 'child_swimming.png' 'couple.png' 'elk.png'...
    'family.png' 'fisher.png' 'foxy.png' 'geckos.png' 'lady_liberty.png' 'lake.png' 'log_seaside.png' 'monument.png' 'native_american.png'...
    'redwood.png' 'roping.png' 'rushmore.png' 'shroom.png' 'snow_leaves.png' 'sunset_sparrow.png' 'sunsetcolor.png' 'swarm.png' 'trolley.png'...
    'turtle.png' 'veggies.png' 'woman.png'};
% dist_num=[5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5];

num_ref=length(img_ref);
num_dist = 116;

images_path_ref = cell(num_ref,1);
images_path_dist = cell(num_dist,1);
to_ref = zeros(num_dist,1);
mos = zeros(num_dist,1);

img_dist=csiq.img_dist;


cnt_dist = 0;
for idx_ref=1:num_ref
    images_path_ref{idx_ref}=[ lib_path, filesep, 'reference', filesep, img_ref{idx_ref} ];
    
%     cum_sum=sum(dist_num(1:idx_ref));
%     sub_start=cum_sum- dist_num(idx_ref)+1;
%     sub_end=cum_sum;
%     for order_dist=sub_start:sub_end;
     for order_dist=1:866;
        name1 = img_dist(order_dist).name_uc;
        str=img_ref{idx_ref};
        len=size(str);
        str1=name1(1:len(2)-5);
        str2=str(1:len(2)-5);
        if(label(order_dist)==6 &&strcmp(str1,str2)==1)
%         name1 = img_dist(order_dist).name_uc;
        name2 = img_dist(order_dist).name;
        
        % find sub directory
        sub_dir = 'contrast';
        
        % find correct path
        dir_path = [ lib_path, filesep, 'distored', filesep, sub_dir ];
        img_path_dist = [ dir_path, filesep, name1 ];
        if ~exist(img_path_dist,'file')
            img_path_dist = [ dir_path, filesep, name2 ];
        end

        cnt_dist = cnt_dist + 1;
        to_ref(cnt_dist) = idx_ref;
        images_path_dist{cnt_dist}=img_path_dist;
        mos(cnt_dist)=dmos(order_dist);
        end
    end
end % for idx_ref

dataset.img_path_ref = images_path_ref;
dataset.img_path_dist = images_path_dist;
dataset.to_ref = to_ref;
dataset.scores = mos;
dataset.scores_type = 'dmos';
end % function getCSIQ



function dataset = getfnoise()
lib_path='dataset/csiq';
csiq=load ('dataset/csiq_dist.mat');
dmos=load( 'dataset/csiq.mat' );
label=dmos.csiq.label;
dmos=dmos.csiq.dmos;


img_ref={ '1600.png' 'aerial_city.png' 'boston.png' 'bridge.png' 'butter_flower.png' 'cactus.png' 'child_swimming.png' 'couple.png' 'elk.png'...
    'family.png' 'fisher.png' 'foxy.png' 'geckos.png' 'lady_liberty.png' 'lake.png' 'log_seaside.png' 'monument.png' 'native_american.png'...
    'redwood.png' 'roping.png' 'rushmore.png' 'shroom.png' 'snow_leaves.png' 'sunset_sparrow.png' 'sunsetcolor.png' 'swarm.png' 'trolley.png'...
    'turtle.png' 'veggies.png' 'woman.png'};
% dist_num=[5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5];

num_ref=length(img_ref);
num_dist = 150;

images_path_ref = cell(num_ref,1);
images_path_dist = cell(num_dist,1);
to_ref = zeros(num_dist,1);
mos = zeros(num_dist,1);

img_dist=csiq.img_dist;


cnt_dist = 0;
for idx_ref=1:num_ref
    images_path_ref{idx_ref}=[ lib_path, filesep, 'reference', filesep, img_ref{idx_ref} ];
    
%     cum_sum=sum(dist_num(1:idx_ref));
%     sub_start=cum_sum- dist_num(idx_ref)+1;
%     sub_end=cum_sum;
%     for order_dist=sub_start:sub_end;
     for order_dist=1:866;
        name1 = img_dist(order_dist).name_uc;
        str=img_ref{idx_ref};
        len=size(str);
        str1=name1(1:len(2)-5);
        str2=str(1:len(2)-5);
        if(label(order_dist)==4 &&strcmp(str1,str2)==1)
%         name1 = img_dist(order_dist).name_uc;
        name2 = img_dist(order_dist).name;
        
        % find sub directory
        sub_dir = 'fnoise';
        
        % find correct path
        dir_path = [ lib_path, filesep, 'distored', filesep, sub_dir ];
        img_path_dist = [ dir_path, filesep, name1 ];
        if ~exist(img_path_dist,'file')
            img_path_dist = [ dir_path, filesep, name2 ];
        end

        cnt_dist = cnt_dist + 1;
        to_ref(cnt_dist) = idx_ref;
        images_path_dist{cnt_dist}=img_path_dist;
        mos(cnt_dist)=dmos(order_dist);
        end
    end
end % for idx_ref

dataset.img_path_ref = images_path_ref;
dataset.img_path_dist = images_path_dist;
dataset.to_ref = to_ref;
dataset.scores = mos;
dataset.scores_type = 'dmos';
end % function getCSIQ


function dataset = getjpeg()
lib_path='dataset/csiq';
csiq=load ('dataset/csiq_dist.mat');
dmos=load( 'dataset/csiq.mat' );
label=dmos.csiq.label;
dmos=dmos.csiq.dmos;


img_ref={ '1600.png' 'aerial_city.png' 'boston.png' 'bridge.png' 'butter_flower.png' 'cactus.png' 'child_swimming.png' 'couple.png' 'elk.png'...
    'family.png' 'fisher.png' 'foxy.png' 'geckos.png' 'lady_liberty.png' 'lake.png' 'log_seaside.png' 'monument.png' 'native_american.png'...
    'redwood.png' 'roping.png' 'rushmore.png' 'shroom.png' 'snow_leaves.png' 'sunset_sparrow.png' 'sunsetcolor.png' 'swarm.png' 'trolley.png'...
    'turtle.png' 'veggies.png' 'woman.png'};
% dist_num=[5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5];

num_ref=length(img_ref);
num_dist = 150;

images_path_ref = cell(num_ref,1);
images_path_dist = cell(num_dist,1);
to_ref = zeros(num_dist,1);
mos = zeros(num_dist,1);

img_dist=csiq.img_dist;


cnt_dist = 0;
for idx_ref=1:num_ref
    images_path_ref{idx_ref}=[ lib_path, filesep, 'reference', filesep, img_ref{idx_ref} ];
    
%     cum_sum=sum(dist_num(1:idx_ref));
%     sub_start=cum_sum- dist_num(idx_ref)+1;
%     sub_end=cum_sum;
%     for order_dist=sub_start:sub_end;
     for order_dist=1:866;
        name1 = img_dist(order_dist).name_uc;
        str=img_ref{idx_ref};
        len=size(str);
        str1=name1(1:len(2)-5);
        str2=str(1:len(2)-5);
        if(label(order_dist)==2 &&strcmp(str1,str2)==1)
%         name1 = img_dist(order_dist).name_uc;
        name2 = img_dist(order_dist).name;
        
        % find sub directory
        sub_dir = 'jpeg';
        
        % find correct path
        dir_path = [ lib_path, filesep, 'distored', filesep, sub_dir ];
        img_path_dist = [ dir_path, filesep, name1 ];
        if ~exist(img_path_dist,'file')
            img_path_dist = [ dir_path, filesep, name2 ];
        end

        cnt_dist = cnt_dist + 1;
        to_ref(cnt_dist) = idx_ref;
        images_path_dist{cnt_dist}=img_path_dist;
        mos(cnt_dist)=dmos(order_dist);
        end
    end
end % for idx_ref

dataset.img_path_ref = images_path_ref;
dataset.img_path_dist = images_path_dist;
dataset.to_ref = to_ref;
dataset.scores = mos;
dataset.scores_type = 'dmos';
end % function getCSIQ



function dataset = getjpeg2000()
lib_path='dataset/csiq';
csiq=load ('dataset/csiq_dist.mat');
dmos=load( 'dataset/csiq.mat' );
label=dmos.csiq.label;
dmos=dmos.csiq.dmos;


img_ref={ '1600.png' 'aerial_city.png' 'boston.png' 'bridge.png' 'butter_flower.png' 'cactus.png' 'child_swimming.png' 'couple.png' 'elk.png'...
    'family.png' 'fisher.png' 'foxy.png' 'geckos.png' 'lady_liberty.png' 'lake.png' 'log_seaside.png' 'monument.png' 'native_american.png'...
    'redwood.png' 'roping.png' 'rushmore.png' 'shroom.png' 'snow_leaves.png' 'sunset_sparrow.png' 'sunsetcolor.png' 'swarm.png' 'trolley.png'...
    'turtle.png' 'veggies.png' 'woman.png'};
% dist_num=[5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5];

num_ref=length(img_ref);
num_dist = 150;

images_path_ref = cell(num_ref,1);
images_path_dist = cell(num_dist,1);
to_ref = zeros(num_dist,1);
mos = zeros(num_dist,1);

img_dist=csiq.img_dist;


cnt_dist = 0;
for idx_ref=1:num_ref
    images_path_ref{idx_ref}=[ lib_path, filesep, 'reference', filesep, img_ref{idx_ref} ];
    
%     cum_sum=sum(dist_num(1:idx_ref));
%     sub_start=cum_sum- dist_num(idx_ref)+1;
%     sub_end=cum_sum;
%     for order_dist=sub_start:sub_end;
     for order_dist=1:866;
        name1 = img_dist(order_dist).name_uc;
        str=img_ref{idx_ref};
        len=size(str);
        str1=name1(1:len(2)-5);
        str2=str(1:len(2)-5);
        if(label(order_dist)==3 &&strcmp(str1,str2)==1)
%         name1 = img_dist(order_dist).name_uc;
        name2 = img_dist(order_dist).name;
        
        % find sub directory
        sub_dir = 'jpeg2000';
        
        % find correct path
        dir_path = [ lib_path, filesep, 'distored', filesep, sub_dir ];
        img_path_dist = [ dir_path, filesep, name1 ];
        if ~exist(img_path_dist,'file')
            img_path_dist = [ dir_path, filesep, name2 ];
        end

        cnt_dist = cnt_dist + 1;
        to_ref(cnt_dist) = idx_ref;
        images_path_dist{cnt_dist}=img_path_dist;
        mos(cnt_dist)=dmos(order_dist);
        end
    end
end % for idx_ref

dataset.img_path_ref = images_path_ref;
dataset.img_path_dist = images_path_dist;
dataset.to_ref = to_ref;
dataset.scores = mos;
dataset.scores_type = 'dmos';
end % function getCSIQ


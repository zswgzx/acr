function lcod(scanner,modality)
% low contrast object detectability test for ACR phantom
% Shengwei Zhang
% Prerequisites: MATLAB R2018b, SPM12, otsu.m
% scanner: uc or mg
% img: name of single axial image for such test
% modality: t1 or t2
% usage: lcod('mg') or lcod('mg','t1'); lcod('uc','t1')
%% check input parameter(s)
if nargin<1
    error('Have to specify scanner site!')
elseif nargin<2
    if contains(scanner,'uc'), error('Have to specify modality for UC scans')
    elseif contains(scanner,'mg')
        img='ACR_T1_AX.nii';modality='t1';
    else, error('No such scanner site')
    end
elseif nargin<3
    if contains(modality,'t1') && contains(scanner,'mg')
        img='ACR_T1_AX.nii';
    elseif contains(modality,'t1') && contains(scanner,'uc')
        img='ACR_Axial_T1.nii';
    elseif contains(modality,'t2') && contains(scanner,'uc')
        img='ACR_Axial_T2_DE_e1.nii';
    else, error('No such combo of scanner + modality')
    end
else, error('Too many input arguments')
end
%% slices extraction
vol=spm_vol(img);
img=spm_read_vols(vol);
img=img(:,:,8:end);
%% find center of spokes
img_seg=otsu(img(:,:,1))-1;
range=7;
bound=zeros(4,1);

% find horizontal coordinates
sum_tmp=sum(img_seg(vol.dim(1)/2-range:vol.dim(1)/2+range,:));
bound(1)=vol.dim(1)/2-find(sum_tmp(vol.dim(1)/2:-1:1)==0,1)+2;
bound(2)=find(sum_tmp(vol.dim(1)/2:end)<range,1)+vol.dim(1)/2-2;

% couldn't detect low contrast objects...
if bound(2)>200, error('Fail to detect low-contrast objects'); end

% in case of otsu failure...
if bound(2)-bound(1)>92 && bound(1)>vol.dim(1)/4, bound(2)=bound(1)+90; 
else, bound(1)=bound(2)-90; 
end
center_col=floor(median(bound(1):bound(2)));

% find vertical coordinates
sum_tmp=img_seg(:,center_col);
bound(3)=vol.dim(1)/2-find(sum_tmp(vol.dim(1)/2:-1:1)==0,1)+2;
bound(4)=find(sum_tmp(vol.dim(1)/2:end)==0,1)+vol.dim(1)/2-2;

% in case of otsu failure...
if bound(4)-bound(3)>92 && bound(4)<vol.dim(1)*.75, bound(3)=bound(4)-90; 
else, bound(4)=bound(3)+90; 
end
center_row=floor(median(bound(3):bound(4)));
%% segment spokes
[x,y]=meshgrid(-center_col+1:vol.dim(1)-center_col,-center_row+1:vol.dim(1)-center_row);
roi=zeros(size(x));
roi((x.^2+y.^2)<diff(bound(1:2))^2/4)=1;
%% save spokes image
img=img.*repmat(roi,[1 1 4]);
new_vol=vol;
new_vol.dim=[size(img,1),size(img,2),size(img,3)];
new_vol.fname=sprintf('lcod-test-%s-spokes.nii',modality);
spm_write_vol(new_vol,img);
%% count spokes
roi=zeros(size(img));
count=zeros(4,1);

for i=4:-1:1
    [count(i),roi(:,:,i)]=count_spoke(img(:,:,i),center_row,center_col,(bound(2)-bound(1))/2,i);
end

% sanity check
for i=1:3
    if count(i+1)<count(i)
        if count(i)==10, count(i+1)=count(i);
        else, count(i+1)=count(i)+1;
        end
    end
end

new_vol.pinfo(1)=1;
new_vol.fname=sprintf('lcod-test-%s-spokes-roi.nii',modality);
spm_write_vol(new_vol,roi);
system('gzip lcod-test*nii');
%% check test criterion
spoke_threshold=37;
if sum(count)>=spoke_threshold
    fprintf('LCOD score = %d, test pass; but double check visually\n',sum(count))
    test_lcod=true;
else
    fprintf('LCOD score = %d, test failed; but double check visually\n',sum(count))
    test_lcod=false;
end
%% write report
fname_result=sprintf('LCOD_results_%s.csv',modality);
fid=fopen(fname_result,'w');
fprintf(fid,'spoke count,result(pass=1;fail=0),slice8 count,slice9 count,slice10 count, slice11 count\n');     
fprintf(fid,'%d,%d,%d,%d,%d,%d\n',sum(count),test_lcod,count(1),count(2),count(3),count(4));
fclose(fid);
end

function [spoke_count,roi]=count_spoke(img,row,col,r,slice)
rads=([1 37.5 74 110 147 183 220 255 289 325]+(slice-1)*9)*pi/180;
spokes=zeros(10,1);
roi=zeros(size(img));

% check each spoke
for i=1:10
    d_row=floor(r*sin(rads(i)));d_col=floor(r*cos(rads(i)));
    
    % https://stackoverflow.com/questions/1940833/how-do-i-create-an-image-matrix-with-a-line-drawn-in-it-in-matlab
    x = [row,row+d_row];
    y = [col,col+d_col];
    nPoints=max(abs(diff(x)),abs(diff(y)))+1;   % Number of points in line
    cIndex=round(linspace(y(1),y(2),nPoints));
    rIndex=round(linspace(x(1),x(2),nPoints));
    index=sub2ind(size(img),rIndex,cIndex);     % Linear indices
    roi(index)=1;
    
    % detect disks
    tmp=img(index);
    tmp=tmp(tmp>0);
    tmp=detrend(tmp);
    tmp(tmp<=0.0001)=0;tmp(tmp>0)=1;
    if length(find_pattern(tmp,[0 1 1 ]))>=2, spokes(i)=1; end
end

% sanity check
tmp=find(spokes,1,'last');
x=find(spokes==0);
if ~isempty(x)
    for i=1:length(x)
        if x(i)<tmp, spokes(x(i))=1; end
    end
end
spoke_count=sum(spokes);
end

function start = find_pattern(array, pattern)
%   Locate a pattern in an array.
%
%   indices = findPattern2(array, pattern) finds the starting indices of
%   pattern within array.
%
%   Example:
%   a = [0 1 4 9 16 4 9];
%   patt = [4 9];
%   indices = findPattern2(a,patt)
%   indices =
%        3     6
%
% https://blogs.mathworks.com/loren/2008/09/08/finding-patterns-in-arrays/

len = length(pattern);
start = find(array==pattern(1));
endVals = start+len-1;
start(endVals>length(array)) = [];

for pattval = 2:len
    % check viable locations in array
    locs = pattern(pattval) == array(start+pattval-1);
    % delete false ones from indices
    start(~locs) = [];
end
end

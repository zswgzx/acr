function spatial_res(scanner,modality)
% spatial resolution test for ACR phantom
% Shengwei Zhang
% Prerequisites: MATLAB R2018b, SPM12, otsu.m, clusterdata.m (build-in)
% scanner: uc or mg
% img: name of single axial image for such test
% modality: t1 or t2
% usage: spatial_res('mg') or spatial_res('mg','t1') ; spatial_res('uc','t1'(or 't2'))
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
%% criteria for such test
% resolutions from left to right 
res=[1.1,1,0.9];
results=false(2,length(res));
%% slice 1 resolution insert extraction
% read volume info
vol=spm_vol(img);
img=spm_read_vols(vol);
img=img(:,:,1);

% assume vol.dim(1)=vol.dim(2)
% roi is empirical...
img_roi=flipud(img(vol.dim(1)/4+1:vol.dim(1)*3/4,...
    vol.dim(1)/4+1:vol.dim(1)*7/16-1)');
img_roi=img_roi(12+1:36-1,32+5:96+9);

img=zeros(23,23,3);
for i=1:3, img(:,:,i)=img_roi(:,23*(i-1)+1:23*i); end
%% save resolution insert
img_roi=img;
img_roi=rot90(img_roi,3);
new_vol=vol;
new_vol.fname=sprintf('spatial-test-%s-res-insert.nii',modality);
new_vol.dim=[23,23,3];
spm_write_vol(new_vol,img_roi);
%% check if resolved
img_roi=otsu(img)-1;
for i=1:3
   % find center hole
   [~, ctr_x]=max(sum(img_roi(:,:,i)));
   [~, ctr_y]=max(sum(img_roi(:,:,i),2));
   
   % check UL & LR
   results(1,i)=resovle_test('ul',img_roi(:,:,i),res(i),ctr_y,ctr_x);
   results(2,i)=resovle_test('lr',img_roi(:,:,i),res(i),ctr_y,ctr_x);
end
%% save binarized resolution insert
new_vol.fname=sprintf('spatial-test-%s-res-insert-bin.nii',modality);
new_vol.pinfo(1)=1;
spm_write_vol(new_vol,rot90(img_roi,3));
system('gzip spatial-test*.nii');
disp('Double check output images visually and edit results as needed!!')
%% write report
fname_result = sprintf('Spatial_resolution_results_%s.csv',modality);
fid = fopen(fname_result,'w');
fprintf(fid,'resolved=1,otherwise=0\n');
fprintf(fid,'1.1mm L-R,1.1mm T-B,1mm L-R,1mm T-B,0.9mm L-R,0.9mm T-B\n');
fprintf(fid,'%d,%d,%d,%d,%d,%d\n', results(:));
fclose(fid);
end

function result=resovle_test(part,img,res,ctr_y,ctr_x)
%% extract ul or lr part
if contains(part,'ul')
    [row,col]=find(img(1:ctr_y,1:ctr_x));
    coord=sortrows([row,col]);
    col=coord(:,2);
    row=coord(:,1);
    c=unique(row);
elseif contains(part,'lr')
    [row,col]=find(img(ctr_y:end,ctr_x:end));
    c=unique(col);
end
%% group c in 4 rows/columns
if length(c)>4
    if contains(part,'ul')
        c=c(histc(row,c)>=4);
    elseif contains(part,'lr')
        c=c(histc(col,c)>=4);
    end
    if length(c)==4, idx=1:4;
    else, idx=clusort(c);
    end
else, idx=1:4;
end
%% check rows/columns
rows=zeros(4,1);
for j=1:4
    idx1=c(idx==j);
    for k=idx1(1):idx1(end)
        if contains(part,'ul')
            if length(find(diff(col(row==k))>=2))>2, rows(j)=1; end
        elseif contains(part,'lr')
            if length(find(diff(row(col==k))>=2))>2, rows(j)=1; end
        end
    end
end
%% show results on screen
for j=1:4
    if contains(part,'ul')
        metric='row';
    elseif contains(part,'lr')
        metric='column';
    end
    if rows(j)==1, fprintf('Insert of %.1fmm: %s %s %d resolved\n',...
            res,upper(part),metric,j);end
end
result=any(rows);
end

function idx=clusort(c)
idx_raw=clusterdata(c,4);
idx=zeros(size(idx_raw));
r=zeros(4,1);
j=1;

% sort cluster indices in ascending order
for i=1:length(idx_raw)
    idx(i)=j;
    if length(find(idx_raw==idx_raw(i)))==1, j=j+1;
    else
        idx_tmp=find(idx_raw==idx_raw(i));
        if i==idx_tmp(end), r(j)=1;j=j+1; end
    end
end

% sanity check for c, watch out for possible improvement
for i=find(r)
    if isempty(i), break;
    idx_tmp=idx==i;
    elseif length(c(idx==i))==2 && diff(c(idx==i))>1
        if i<4, idx(idx_tmp(end))=i+1;
        else, idx(idx_tmp(1))=i-1;
        end
    elseif length(c(idx==i))==3
        if i>1 && ((idx_tmp(2)+2<=size(c,1) && c(idx_tmp(2)+2)-c(idx_tmp(2))>1) || idx_tmp(2)+2>size(c,1))
            idx(idx_tmp(1))=i-1;
        end
    elseif length(c(idx==i))>=4
        if isempty(find(diff(c(idx==i))>1, 1)) && i>1
            idx(idx_tmp(1))=i-1;idx(idx_tmp(2))=i-1;
        end
    end
end
end

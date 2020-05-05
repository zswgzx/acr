function slice_pos(scanner,modality)
% slice position test for ACR phantom
% Shengwei Zhang
% Prerequisites: MATLAB R2018b, SPM12, otsu.m
% scanner: uc or mg
% img: name of single axial image for such test
% modality: t1 or t2
% usage: slice_pos('mg') or slice_pos('mg','t1'); slice_pos('uc','t1')
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
%% read volume/header info
vol=spm_vol(img);
hdr=spm_read_hdr(img);
len=hdr.dime.pixdim(2);
%% slice 1 & 11 bar extraction
img=spm_read_vols(vol);
img(:,:,2)=img(:,:,11);img=img(:,:,1:2);
bar_range=10;
bar_diff=zeros(2,1);

% aassume t1_vol.dim(1)=t1_vol.dim(2) & t1_hdr.dime.pixdim(2)=t1_hdr.dime.pixdim(3)
for i=1:2
    if i==1, bar_diff(i)=run(vol.dim(1),img(:,:,i),bar_range,vol,modality,i,len);
    else, bar_diff(i)=run(vol.dim(1),img(:,:,i),bar_range,vol,modality,11,len);
    end
end
%% check test criterion
bar_len_limit=[4 7];
if max(abs(bar_diff))<=bar_len_limit(2)
    fprintf('Max. abs. bar len. difference is %.2fmm, test passed\n',max(abs(bar_diff)))
    test_stat = true;
else
    fprintf('Max. abs. bar len. difference is %.2fmm, test failed\n',max(abs(bar_diff)))
    test_stat = false;
end
%% write report
fname_result=sprintf('Slice_position_results_%s.csv',modality);
fid= fopen(fname_result,'w');
fprintf(fid,'Max. abs. bar len. difference (mm),result(pass=1;fail=0),bar len. diff. slice1 (mm),bar len. diff. slice11 (mm)\n');
fprintf(fid,'%.2f,%d,%.2f,%.2f\n',max(abs(bar_diff)),test_stat,bar_diff(1),bar_diff(2));
fclose(fid);
end

function bar_diff=run(dim,img,range,vol,modality,slice,len)
bound=zeros(2,1);
%% rough bar segmentation
img=flipud(img(dim/2-range:dim/2+range,(0.75*dim-range):end)');
img_seg=otsu(img)-1;
%% improve segmentation
sum_ver=sum(img_seg,2);
bound(1)=find(sum_ver,1);
tmp=find(sum_ver==2*range+1);
bound(2)=tmp(2);
if bound(2)+3<=size(img,1)
    img=img(bound(1):bound(2)+3,:);
else, img=img(bound(1):end,:);
end
img_seg=otsu(img)-1;
%% save bar image
new_vol=vol;
new_vol.dim=[size(img,2),size(img,1),1];
new_vol.fname=sprintf('slice-position-test-%s-s%d.nii',modality,slice);
spm_write_vol(new_vol,rot90(img,3));

new_vol.pinfo(1)=1;
new_vol.fname=sprintf('slice-position-test-%s-s%d-bin.nii',modality,slice);
spm_write_vol(new_vol,rot90(img_seg,3));
system(sprintf('gzip slice-position-test-%s-s%d*nii',modality,slice));
%% calculate bar length
sum_hor=sum(img_seg(floor(size(img_seg,1)/2):end,:));
sum_hor=sum_hor(sum_hor<1.5*range);
[counts,lens]=hist(sum_hor,unique(sum_hor));
lens=lens(counts>1);
counts=counts(counts>1);
if length(lens)==4
    [counts_ord,idx]=sort(counts,'descend');
    if counts_ord(1)>counts_ord(2)
        bound(1)=lens(idx(1));
        if counts_ord(2)>counts_ord(3)
            bound(2)=lens(idx(2));
        elseif counts_ord(3)>counts_ord(4)
            bound(2)=lens(idx(3));
        else, bound(2)=lens(idx(4));
        end
    end
    bar_diff=-len*diff(bound);
elseif length(lens)==3
    if lens(2)-lens(1)<2
        bound(2)=lens(3);
        if counts(1)>=counts(2), bound(1)=lens(1);
        else, bound(1)=lens(2);
        end
    else
        bound(1)=lens(1);
        if counts(2)<=counts(3), bound(2)=lens(3);
        else, bound(2)=lens(2);
        end
    end
    bar_diff=-len*diff(bound);
elseif length(lens)==2
    bar_diff=-len*diff(lens);
elseif length(lens)==1
    bar_diff=-len;
end
end

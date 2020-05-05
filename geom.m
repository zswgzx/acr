function geom(scanner)
% geometry test for ACR phantom
% Shengwei Zhang, RADC
% Prerequisites: MATLAB R2018b, SPM12, otsu.m
% -scanner: mg or uc
% -loc: name of single localizer image for such test
% -t1: name of single t1 axial image for such test
% usage: geom('mg') or geom('uc')
%% check input parameter(s)
if nargin<1
    error('Have to specify scanner site!')
elseif nargin<2
    if contains(scanner,'uc')
        loc='ACR_Sagittal_Locator.nii';
        t1='ACR_Axial_T1.nii';
    elseif contains(scanner,'mg')
        loc='ACR_T1_SAG_h-f.nii';
        t1='ACR_T1_AX.nii';
    else, error('No such scanner site')
    end
else, error('Too many input arguments')
end
%% criteria for such test
thr_range=3;    % threshold of acceptable range in mm
loc_thr=148;    % inside end-to-end length of phantom in mm
t1_thr=190;     % inside diameter of phantom in mm
%% Localizer test
% read volume/header info
vol=spm_vol(loc);
hdr=spm_read_hdr(loc);

% binarize image from original
img=spm_read_vols(vol);
seg=otsu(img)-1;
% figure;imshow(seg,[]);

% zero margins due to possible artifacts
margin=10;
seg(1:margin,:)=0;
seg(vol.dim(1)-margin:vol.dim(1),:)=0;
seg(:,1:margin)=0;
seg(:,vol.dim(2)-margin:vol.dim(2))=0;

% sum along rows in mm
sum_tmp=sum(seg,2)*hdr.dime.pixdim(3);
loc_range=sum_tmp(sum_tmp >=loc_thr-thr_range);

if ~isempty(loc_range) && abs(max(loc_range)-loc_thr)<=thr_range
    disp('Localizer passed the test');
    test_loc_hor=true;
elseif isempty(loc_range)
    % in case of otsu failure...
    thr=ceil(mean(img(seg==0))+std(img(seg==0)));
    seg=img;
    seg(seg<=thr)=0;
    seg(seg>thr)=1;
    
    seg(1:margin,:)=0;
    seg(vol.dim(1)-margin:vol.dim(1),:)=0;
    seg(:,1:margin)=0;
    seg(:,vol.dim(2)-margin:vol.dim(2))=0;
    
    sum_tmp=sum(seg,2)*hdr.dime.pixdim(3);
    loc_range=sum_tmp(sum_tmp >=loc_thr-thr_range);
    if abs(max(loc_range)-loc_thr)<=thr_range
        disp('Localizer passed the test');
        test_loc_hor=true;
    else
        disp('Localizer failed the test');
        test_loc_hor=false;
    end
else
    disp('Localizer failed the test');
    test_loc_hor=false;
end

% save binarized image
bin_vol=vol;
bin_vol.fname='geom-test-localizer-bin.nii';
spm_write_vol(bin_vol,seg);
%% T1 test
% read volume/header info
vol=spm_vol(t1);
hdr=spm_read_hdr(t1);

% binarize slice 1 & 5
img=spm_read_vols(vol);
img(:,:,2)=img(:,:,5);img=img(:,:,1:2);
seg=zeros(size(img));
test_hor=zeros(2,1);test_ver=test_hor;
disc_hor=test_hor;disc_ver=test_hor;

for i=1:2
    img_tmp=img(:,:,i);
    seg_tmp=otsu(img_tmp)-1;

    % zero margins due to possible artifacts
    seg_tmp(1:margin,:)=0;
    seg_tmp(vol.dim(1)-margin:vol.dim(1),:)=0;
    seg_tmp(:,1:margin)=0;
    seg_tmp(:,vol.dim(2)-margin:vol.dim(2))=0;
    seg(:,:,i)=seg_tmp;

    % run test on horizontal & vertical directions
    disc_hor(i)=hor_ver_test(seg_tmp,hdr.dime.pixdim(2),'hor');
    disc_ver(i)=hor_ver_test(seg_tmp,hdr.dime.pixdim(3),'ver');
    
    % in case otsu fails...
    if disc_ver(i)>t1_thr+thr_range || disc_ver(i)<t1_thr-thr_range
        thr=ceil(mean(img_tmp(seg_tmp==0))+std(img_tmp(seg_tmp==0)));
        seg_tmp=img_tmp;
        seg_tmp(seg_tmp<=thr)=0;
        seg_tmp(seg_tmp>thr)=1;
        
        seg_tmp(1:margin,:)=0;
        seg_tmp(vol.dim(1)-margin:vol.dim(1),:)=0;
        seg_tmp(:,1:margin)=0;
        seg_tmp(:,vol.dim(2)-margin:vol.dim(2))=0;
        disc_ver(i)=hor_ver_test(seg_tmp,hdr.dime.pixdim(2),'ver');
        
        if disc_ver(i)<=t1_thr+thr_range && disc_ver(i)>=t1_thr-thr_range
            seg(:,:,i)=seg_tmp;
            disc_hor(i)=hor_ver_test(seg_tmp,hdr.dime.pixdim(2),'hor');
        end
    end
    
    if i==1
        test_hor(i)=len_test(i,disc_hor(i),t1_thr,thr_range,'horizontal');
        test_ver(i)=len_test(i,disc_ver(i),t1_thr,thr_range,'vertical');
    else
        test_hor(i)=len_test(5,disc_hor(i),t1_thr,thr_range,'horizontal');
        test_ver(i)=len_test(5,disc_ver(i),t1_thr,thr_range,'vertical');
    end
    
    % save binarized images
    bin_vol=vol;
    if i==1, bin_vol.fname=sprintf('geom-test-slice%d-bin.nii',i);
    else, bin_vol.fname=sprintf('geom-test-slice%d-bin.nii',5);
    end
    bin_vol.dim(3)=1;
    spm_write_vol(bin_vol,seg(:,:,i));
end
system('gzip geom-*.nii');

% run test on 45 degree left & right directions
seg=seg(:,:,2);
idx=zeros(vol.dim(1),2);
test_diag=zeros(2,1);disc_diag=test_diag;
for i=1:2
    for j=1:vol.dim(1)
        if i==1, margin=find(diag(seg,j-vol.dim(1)/2-1),1,'first');
        else, margin=find(diag(fliplr(seg),j-vol.dim(1)/2-1),1,'first');
        end

        if isempty(margin), continue
        else
            idx(j,1)=margin;
            if i==1, idx(j,2)=find(diag(seg,j-vol.dim(1)/2-1),1,'last');
            else, idx(j,2)=find(diag(fliplr(seg),j-vol.dim(1)/2-1),1,'last');
            end
        end
    end
    disc_diag(i)=max(diff(idx,1,2)+1)*hdr.dime.pixdim(3)*sqrt(2);
    
    if i==1
        test_diag(i)=len_test(5,disc_diag(i),t1_thr,thr_range,'45 deg. left');
    else
        test_diag(i)=len_test(5,disc_diag(i),t1_thr,thr_range,'45 deg. right');
    end
end
%% write report
fname_result = sprintf('Geometry_results.csv');
fid = fopen(fname_result,'w');
fprintf(fid,'Parameter,Max. Distance (mm),Result(pass=1;fail=0)\n');     
fprintf(fid,'localizer               ,%.2f,%d\n',max(loc_range),test_loc_hor);
fprintf(fid,'T1 slice 1 horizontal   ,%.2f,%d\n',disc_hor(1),test_hor(1));
fprintf(fid,'T1 slice 1 vertical     ,%.2f,%d\n',disc_ver(1),test_ver(1));
fprintf(fid,'T1 slice 5 horizontal   ,%.2f,%d\n',disc_hor(2),test_hor(2));
fprintf(fid,'T1 slice 5 vertical     ,%.2f,%d\n',disc_ver(2),test_ver(2));
fprintf(fid,'T1 slice 5 45 deg. left ,%.2f,%d\n',disc_diag(1),test_diag(1));
fprintf(fid,'T1 slice 5 45 deg. right,%.2f,%d\n',disc_diag(2),test_diag(2));
fclose(fid);
end

function res=len_test(slice,max,thr,range,axis)
if abs(max-thr) <= range
    fprintf('T1 measurement slice %d passed the test (%s line)\n',slice,axis)
    res=true;
else
    fprintf('T1 measurement slice %d failed the test (%s line)\n',slice,axis)
    res=false;
end
end

function len=hor_ver_test(img,length,dir)
len=zeros(size(img,1),2);
for i = 1:size(img,1)
    if contains(dir,'ver')
        if isempty(find(img(i,:),1,'first'))
            continue
        else
            len(i,1)=find(img(i,:),1,'first');
            len(i,2)=find(img(i,:),1,'last');
        end
    elseif contains(dir,'hor')
        if isempty(find(img(:,i),1,'first'))
            continue
        else
            len(i,1)=find(img(:,i),1,'first');
            len(i,2)=find(img(:,i),1,'last');
        end
    else, error('Either hor(izontal) or ver(tical) should be selected')
    end
end
len=max(diff(len,1,2)+1)*length;
end

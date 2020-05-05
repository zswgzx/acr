function slice_th(scanner,modality)
% slice thickness test for ACR phantom
% Shengwei Zhang
% Prerequisites: MATLAB R2018b, SPM12, otsu.m
% scanner: uc or mg
% img: name of single axial image for such test
% modality: t1 or t2
% usage: slice_th('mg') or slice_th('mg','t1'); slice_th('uc','t1'(or 't2'))
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
thickness_std=5;
range_std=[0.7 1];
%% read volume/header info
vol=spm_vol(img);
hdr=spm_read_hdr(img);
len=hdr.dime.pixdim(2);
%% slice 1 ramps extraction
img=spm_read_vols(vol);
img=img(:,:,1);
ramp_range=10;

% assume t1_vol.dim(1)=t1_vol.dim(2) & t1_hdr.dime.pixdim(2)=t1_hdr.dime.pixdim(3)
% rough ramp extraction
img_roi=img(:,vol.dim(1)/2-ramp_range:vol.dim(1)/2+ramp_range)';
img_seg=otsu(img_roi)-1;

% improve ramp extraction
sum_tmp=sum(img_seg);
bound=find(sum_tmp==2*ramp_range+1,1);
img_seg(:,1:bound)=1;
bound=find(sum_tmp==2*ramp_range+1,1,'last');
img_seg(:,bound:end)=1;
img_roi(img_seg==1)=0;

threshold=mean(img_roi(img_roi>0))+2*std(img_roi(img_roi>0));
img_roi(img_roi>=threshold)=0;
img_seg=otsu(img_roi)-1;
%% save ramps image
new_vol=vol;
new_vol.dim=[size(img,1),size(img,2),1];
new_vol.fname=sprintf('slice-thickness-test-%s.nii',modality);
spm_write_vol(new_vol,rot90(img_roi,3));
system(sprintf('gzip %s',new_vol.fname));

new_vol.pinfo(1)=1;
new_vol.fname=sprintf('slice-thickness-test-%s-bin.nii',modality);
spm_write_vol(new_vol,rot90(img_seg,3));
system(sprintf('gzip %s',new_t1_vol.fname));
%% main process
lens=zeros(2,1);
tmp=lens;
if contains(scanner,'mg')
    % separate top down bar
    sum_tmp=sum(img_seg,2);
    bound=[find(sum_tmp,1)+2,find(sum_tmp,1,'last')-2];
    mid=find(sum_tmp==min(sum_tmp(bound(1):bound(2))));
    
    % calculate lengths for each bar
    for i=1:2
        if i==1
            roi=img_seg(find(sum_tmp,1):mid-1,:);
        else
            roi=img_seg(mid+1:find(sum_tmp,1,'last'),:);
        end
        
        % find boundary of left/right
        sum_tmp=fliplr(sum(roi(:,1:t1_vol.dim(1)/2)));
        tmp(1)=t1_vol.dim(1)/2+2-find(sum_tmp==1,1);
        tmp(2)=t1_vol.dim(1)/2+1-find(sum_tmp==2,1);
        bound(1)=ceil(median(tmp(1):tmp(2)));
        
        sum_tmp=sum(roi(:,t1_vol.dim(1)/2+1:end));
        tmp(1)=t1_vol.dim(1)/2+find(sum_tmp==2,1);
        tmp(2)=t1_vol.dim(1)/2+find(sum_tmp==1,1);
        bound(2)=floor(median(tmp(1):tmp(2)));
        lens(i)=diff(bound)+1;
        
        if lens(i)>60
            % check bound(1)
            
            % check bound(2)
            tmp(1)=t1_vol.dim(1)/2+find(sum_tmp==3,1);
            tmp(2)=bound(2);
            bound(2)=floor(median(tmp(1):tmp(2)));
        end
        lens(i)=diff(bound)+1;
    end
    sl_thickness=0.2*len*(lens(1)*lens(2))/(lens(1)+lens(2));
else
    % can't separate top/down bar, just find boundary of left/right
    sum_tmp=fliplr(sum(img_seg(:,1:t1_vol.dim(1)/2)));
    threshold=sort(unique(sum_tmp),'descend');
    bound(1)=t1_vol.dim(1)/2+2-find(sum_tmp==threshold(3),1);
    
    sum_tmp=sum(img_seg(:,t1_vol.dim(1)/2+1:end));
    threshold=sort(unique(sum_tmp),'descend');
    bound(2)=t1_vol.dim(1)/2+find(sum_tmp==threshold(3),1);
    lens(1)=diff(bound)+1;
    lens(2)=lens(1);
    sl_thickness=0.1*len*lens(1);
end
%% check slice thickness
fprintf('top = %.2f mm, bottom = %.2f mm\n',lens(1)*len,lens(2)*len)

if abs(thickness_std-sl_thickness)<=range_std(1)
    fprintf('Slice thickness (%.2fmm) test pass\n',sl_thickness)
    test_stat=true;
elseif abs(thickness_std-sl_thickness)>=range_std(2)
    fprintf('Slice thickness (%.2fmm) test fail\n',sl_thickness)
    test_stat=false;
else
    fprintf('Slice thickness (%.2fmm) test ambiguous, check outputs\n',sl_thickness)
    test_stat=.5;
end
disp('Double check output images visually and edit results as needed!!')
%% write report
fname_result = sprintf('Slice_thickness_results_%s.csv',modality);
fid= fopen(fname_result,'w');
fprintf(fid,'top (mm),bottom (mm),slice thickness (mm),result(pass=1;fail=0;unclear=.5)\n');
fprintf(fid,'%.2f,%.2f,%.2f,%.1f\n',lens(1)*len,lens(2)*len,sl_thickness,test_stat);
fclose(fid);

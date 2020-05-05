function iiu_psg(scanner,modality)
% image intensity uniformity & percent. signal ghosting test for ACR phantom
% Shengwei Zhang
% Prerequisites: MATLAB R2018b, SPM12, otsu.m
% scanner: uc or mg
% img: name of single axial image for such test
% modality: t1 or t2
% usage: iiu_psg('mg') or iiu_psg('mg','t1'); iiu_psg('uc','t1')
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
%% slice 7 extraction
img=spm_read_vols(vol);
img=img(:,:,7);
img_seg=otsu(img)-1;
%% detect dark square
range=7;
sum_tmp=sum(img_seg(vol.dim(1)/2-range:vol.dim(1)/2+range,:));
sum_tmp(sum_tmp==2*range+1)=0;
bound_rad=find(sum_tmp,1)-1;
% avoid false positive
if bound_rad<vol.dim(1)/2
    bound_rad=find(sum_tmp,2);
    bound_rad=bound_rad(end)-1; 
end
%% find center and 4 rectangular rois of phantom
width=16;
sum_tmp=sum(img_seg);
bound=zeros(1,4);
for i=1:2
    if i==1, bound(i)=find(sum_tmp,1);
    else, bound(i)=find(sum_tmp,1,'last');
    end
end
center_col=floor(median(bound(1):bound(2)));

sum_tmp=sum(img_seg,2);
for i=1:2
    if i==1, bound(i+2)=find(sum_tmp,1);
    else, bound(i+2)=find(sum_tmp,1,'last');
    end
end
center_row=floor(median(bound(1):bound(2)));

img_rois=zeros(size(img,1),size(img,2),4);
if vol.dim(1)-bound(2)<width
    fprintf('Top ROI would probably be different from the others, double check results\n')
    img_rois(center_row-4*width:center_row+4*width,bound(2)+2:vol.dim(1)-1,1)=1;
else
    range=floor(median(bound(2)+2:vol.dim(1)-1));
    img_rois(center_row-4*width:center_row+4*width,range-width/2+1:range+width/2,1)=1;
end

if bound(1)<width
    fprintf('Bottom ROI would probably be different from the others, double check results\n')
    img_rois(center_row-4*width:center_row+4*width,2:bound(1)-2,2)=1;
else
    range=ceil(median(2:bound(1)-2));
    img_rois(center_row-4*width:center_row+4*width,range-width/2+1:range+width/2,2)=1;
end

if vol.dim(1)-bound(4)<width
    fprintf('Right ROI would probably be different from the others, double check results\n')
    img_rois(bound(4)+2:vol.dim(1)-1,center_col-4*width:center_col+4*width,3)=1;
else
    range=floor(median(bound(4)+2:vol.dim(1)-1));
    img_rois(range-width/2+1:range+width/2,center_col-4*width:center_col+4*width,3)=1;
end

if bound(3)<width
    fprintf('Left ROI would probably be different from the others, double check results\n')
    img_rois(2:bound(3)-2,center_col-4*width:center_col+4*width,4)=1;
else
    range=ceil(median(2:bound(3)-2));
    img_rois(range-width/2+1:range+width/2,center_col-4*width:center_col+4*width,4)=1;
end
%% define large roi boundary
[x,y]=meshgrid(-center_col+1:vol.dim(1)-center_col,-center_row+1:vol.dim(1)-center_row);
roi_large=zeros(size(x));
roi_large((x.^2+y.^2)<(bound_rad-center_col-1)^2)=1;
img_roi=img.*roi_large;

if sum(roi_large(:))*len^2>2.05*10^4 || sum(roi_large(:))*len^2<1.95*10^4
   fprintf('Area of chosen large ROI is %.2f mm^2\n',sum(roi_large(:))*len^2) 
end
%% calculate ghosting measurement
mean_large=mean(img(roi_large>0));
mean_top=mean(img(img_rois(:,:,1)>0));
mean_bottom=mean(img(img_rois(:,:,2)>0));
mean_right=mean(img(img_rois(:,:,3)>0));
mean_left=mean(img(img_rois(:,:,4)>0));
%% save large roi image
new_vol=vol;
new_vol.dim=[size(img,1),size(img,2),1];
new_vol.fname=sprintf('intensity-uniform-test-%s-roiL.nii',modality);
spm_write_vol(new_vol,img_roi);
%% define small rois boundary
tmp=unique(img_roi(:));tmp=tmp(2:end);

% min of tentative 'low-signal' intensity
x=max(tmp(tmp<floor(mean(img_roi(img_roi>0))-2*std(img_roi(img_roi>0)))));
if ~isempty(x)
    x=tmp(tmp<floor(mean(img_roi(img_roi>0))-2*std(img_roi(img_roi>0))));
    while length(find(img_roi==x(1)))>1
        x=x(2:end);
        if length(x)==1, break; end
    end
    x=find(img_roi==x(1));
else, x=find(img_roi==tmp(1));
end
if length(x)>1, y=randperm(length(x));x=x(y(1)); end

% max of tentative 'high-signal' intensity
y=tmp(tmp>floor(mean(img_roi(img_roi>0))+2*std(img_roi(img_roi>0))));
if ~isempty(y)
    while length(find(img_roi==tmp(end)))>1
        y=y(1:end-1);
        if length(y)==1, break; end
    end
end
y=find(img_roi==tmp(end));
if length(y)>1, tmp=randperm(length(y));y=y(tmp(1)); end

% find coordinates of such intensities
tmp=zeros(2,2);
tmp(:,1)=[ceil(x/vol.dim(1));rem(x,vol.dim(1))];
tmp(:,2)=[ceil(y/vol.dim(1));rem(y,vol.dim(1))];

% find small rois; note that if the max/min is at the edge of large roi,
% small roi of that region may not be fully in large roi but it's mostly fine
img_roi_sm=zeros(size(img,1),size(img,2),2);
len=5;
for i=1:2
    [x,y]=meshgrid(-tmp(1,i)+1:vol.dim(1)-tmp(1,i),-tmp(2,i)+1:vol.dim(1)-tmp(2,i));
    bound=zeros(size(x));
    bound((x.^2+y.^2)<=(len+1)^2)=1;
    x=bound+roi_large;
    x(x>0)=1;
    if  sum(x(:))>sum(roi_large(:))
        if tmp(1,i)>vol.dim(1)/2 && tmp(2,i)<vol.dim(1)/2
            tmp(1,i)=tmp(1,i)-len+1;tmp(2,i)=tmp(2,i)+len-1;
        elseif tmp(1,i)<vol.dim(1)/2 && tmp(2,i)<vol.dim(1)/2
            tmp(1,i)=tmp(1,i)+len-1;tmp(2,i)=tmp(2,i)+len-1;
        elseif tmp(1,i)<vol.dim(1)/2 && tmp(2,i)>vol.dim(1)/2
            tmp(1,i)=tmp(1,i)+len-1;tmp(2,i)=tmp(2,i)-len+1;
        else, tmp(1,i)=tmp(1,i)-len+1;tmp(2,i)=tmp(2,i)-len+1;
        end
        [x,y]=meshgrid(-tmp(1,i)+1:vol.dim(1)-tmp(1,i),-tmp(2,i)+1:vol.dim(1)-tmp(2,i));
        bound=zeros(size(x));
        bound((x.^2+y.^2)<=(len+1)^2)=1;
    end
    img_roi_sm(:,:,i)=img.*bound;
    sum_tmp=img_roi_sm(:,:,i);
    if i==1, mean_low=mean(sum_tmp(sum_tmp>0));
    else, mean_high=mean(sum_tmp(sum_tmp>0));
    end
end
%% save small roi images
new_vol.fname=sprintf('intensity-uniform-test-%s-roiSl.nii',modality);
spm_write_vol(new_vol,img_roi_sm(:,:,1));
new_vol.fname=sprintf('intensity-uniform-test-%s-roiSh.nii',modality);
spm_write_vol(new_vol,img_roi_sm(:,:,2));
system('gzip intensity-uniform-test*nii');
%% check test criterion
% intensity uniformity
piu_3t_threshold=[.8 .82];
piu=1-(mean_high-mean_low)/(mean_high+mean_low);
if piu>=piu_3t_threshold(2)
    fprintf('PIU = %.2f, test pass\n',piu)
    test_piu=true;
elseif piu<=piu_3t_threshold(1)
    fprintf('PIU = %.2f, test failed\n',piu)
    test_piu=false;
else
    fprintf('PIU = %.2f, test unclear\n',piu)
    test_piu=.5;
end

% ghosting ratio
ratio_threshold=.025;
ghost_ratio = abs(((mean_top+mean_bottom)-(mean_left+mean_right)))/(2*mean_large);
if ghost_ratio<=ratio_threshold
    fprintf('Ghost ratio is %.4f, test pass\n',ghost_ratio)
    test_psg=true;
else
    fprintf('Ghost ratio is %.4f, test failed\n',ghost_ratio)
    test_psg=false;
end
%% write report
fname_result=sprintf('PIU_results_%s.csv',modality);
fid=fopen(fname_result,'w');
fprintf(fid,'PIU,result(pass=1;fail=0;unclear=.5),high,low\n');     
fprintf(fid,'%.2f,%.1f,%.2f,%.2f\n',piu,test_piu,mean_high,mean_low);
fclose(fid);

fname_result=sprintf('PSG_results_%s.csv',modality);
fid=fopen(fname_result,'w');
fprintf(fid,'PSG,result(pass=1;fail=0),top,bottom,right,left\n');     
fprintf(fid,'%.4f,%d,%.2f,%.2f,%.2f,%.2f\n',ghost_ratio,test_psg,mean_top,mean_bottom,...
    mean_right,mean_left);
fclose(fid);

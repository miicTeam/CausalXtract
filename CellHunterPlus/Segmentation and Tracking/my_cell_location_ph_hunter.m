function [data_sp,data_tu]=my_cell_location_ph_hunter(vidFrames,r_min_sp)

close all
fac=1;
m = 1;
Nframes = size(vidFrames,3);
data_sp=[];
data_tu = [];
radii_sp = [];
polarity='bright';

%r_min_sp=4;
r_min_tu_int=14/m;

for i = 1 : Nframes
    [i Nframes]

    Ib2_sp=im2double(vidFrames(:,:,i));
    Ib2_sp(Ib2_sp>1)=1;
    Ib2_sp(Ib2_sp<0)=0;
    % tracking tumor cell (one cell into the center of the frame)
    [~ ,c_tu_int,r_tu_int]=my_CHT_cell(Ib2_sp,polarity,round(max(6/m,r_min_tu_int*fac-r_min_tu_int*fac*0.10)),round(max(14/m,r_min_tu_int*fac+r_min_tu_int*fac*0.10)),fac);
    data_tu = [data_tu; repmat(i,length(r_tu_int),1) c_tu_int(:,[1 , 2])]; %one tu cell
    c_tu=[c_tu_int];
    r_tu = [r_tu_int];
    x = 1:size(Ib2_sp,1);
    y = 1:size(Ib2_sp,2);
    [Y,X] = meshgrid(y,x);
    mask = Ib2_sp;
    for s = 1 : size(c_tu,1)
        mask((X-c_tu(s,2)).^2+(Y-c_tu(s,1)).^2 <= (r_tu(s)).^2) = 0;

    end
    Ib3_sp= mask;
    [~ ,c_sp,r_sp]=my_CHT_cell(Ib3_sp,polarity,round(max(1,r_min_sp*fac-r_min_sp*fac*0.15)),round(max(1,r_min_sp*fac+r_min_sp*fac*0.15)),fac);

    if exist('c_sp') && not(isempty(c_sp))
        data_sp=[data_sp;repmat(i,length(r_sp),1) c_sp(:,[1 , 2])];
    end


    %
    if mod(i,500)==0
        figure(1);
        imshow(Ib3_sp); hold on;
        viscircles(c_sp,r_sp,'EdgeColor','r');
        %
        % %hold on;
        % %viscircles(c_tu,r_tu,'EdgeColor','g');
        drawnow;
        hold off;
        pause(0.5);
    end


    clear c_sp r_sp c_tu r_tu Ib Ib2_sp Ib2_tu Ib3_sp CHT c_sp2 r_sp2 mask c_tu_max c_tu_min c_tu_int r_tu_max r_tu_min r_tu_int;
end


id_sp=data_sp(1,1);
if id_sp > 1
    id_d=data_sp(data_sp(:,1)==id_sp,:);
    for i = 1 : id_sp-1
        data_sp2((i-1)*size(id_d,1)+1:i*size(id_d,1),2:3)=data_sp(1:size(id_d,1),2:3);
        data_sp2((i-1)*size(id_d,1)+1:i*size(id_d,1),1)=i;
    end
    data_sp2(size(data_sp2,1)+1:size(data_sp2,1)+size(data_sp,1),:)=data_sp;
    clear data_sp
    data_sp=data_sp2;
    clear data_sp2 id_sp id_d
end
end
function [CHT,centers,radii]=my_CHT_cell(B,polarity,rmin,rmax,fac)
warning off;
% fac=1;
if fac>1
I2=imresize(B,fac,'bilinear');
%I2=my_morph_sharp(I2,round(rmin/2));
else
I2=B;
end
centers =[];
sens=0.96;%0.95;
[centers,radii,metric]=imfindcircles(I2,[rmin rmax],'ObjectPolarity',polarity,'Sensitivity',sens);
while isempty(centers) && sens<0.99
sens=min(0.99,sens+0.01);
[centers,radii,metric]=imfindcircles(I2,[rmin rmax],'ObjectPolarity',polarity,'Sensitivity',sens);
end
alpha=linspace(0,360,361)/180*pi;
CHT=zeros(size(I2));
for k = 1 : length(radii)
    x=min(size(I2,1),max(1,floor(centers(k,2)+radii(k)*sin(alpha))));
    y=min(size(I2,2),max(1,floor(centers(k,1)+radii(k)*cos(alpha))));
    idx=sub2ind(size(I2),x,y);
    CHT(idx)=1;
    clear x y;
end
CHT=imfill(CHT,'holes');
CHT=imresize(CHT,1/fac,'nearest');
centers=centers/fac;
radii=radii/fac;
end







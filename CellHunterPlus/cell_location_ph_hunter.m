function [data_sp,data_tu]=cell_location_ph_hunter(vidFrames,...
    r_min_tu_int, r_min_sp, polarity, dist_tu, flag_imm)

close all
fac=1;
m=1;
Nframes =size(vidFrames,3);
r_min_sp=r_min_sp/m;
r_min_tu_int=r_min_tu_int/m;

data_sp=[];
data_tu=[];
x_c=size(vidFrames,1)/2;
y_c=size(vidFrames,2)/2;

for i = 1 : Nframes
    [i Nframes]
    Ib2_sp=im2double(vidFrames(:,:,i));
    Ib2_sp(Ib2_sp>1)=1;
    Ib2_sp(Ib2_sp<0)=0;
    %Ib2_sp=imadjust(Ib2_sp,stretchlim(Ib2_sp));
    roi=imcrop(Ib2_sp,[y_c-50 x_c-50 100 100]);

    % tracking tumor cell (one cell in the center of the frame)
    % localize cancer cells in the whole frame
    [~ ,c_tu_int,r_tu_int]=my_CHT_cell(Ib2_sp,polarity,round(max(6/m,r_min_tu_int*fac-r_min_tu_int*fac*0.10)),round(min(16/m,r_min_tu_int*fac+r_min_tu_int*fac*0.10)),fac,0.99);
    % localize on cancer cell in the middle of the frame
    [~ ,c_tu_int_c,r_tu_int_c]=my_CHT_cell(roi,polarity,round(max(6/m,r_min_tu_int*fac-r_min_tu_int*fac*0.10)),round(min(16/m,r_min_tu_int*fac+r_min_tu_int*fac*0.10)),fac,0.99);
    % if more than two cells are identified in the middle of the frame,
    % for i=1 (first frame), keep the cell closer to the center,
    % as long as the distance from the center is less than dist_tu px;
    % for i>1, keep the cell closer to the one previously identified in
    % frame i-1.
    if exist('c_tu_int_c', 'var') && not(isempty(c_tu_int_c))
        n_cells=size(c_tu_int_c,1);
        if n_cells>=1
            dist_c=zeros(n_cells,1);
            for num_c=1:n_cells
                dist_c(num_c)=sqrt((c_tu_int_c(num_c,2)-size(roi,1)/2).^2+(c_tu_int_c(num_c,1)-size(roi,2)/2).^2);
            end
            final_idx=dist_c<dist_tu;
            dist_c=dist_c(final_idx,:);
            c_tu_int_c=c_tu_int_c(final_idx,:);
            r_tu_int_c=r_tu_int_c(final_idx);
            if size(c_tu_int_c,1)>1
                if i==1
                    [~,I]=min(dist_c);
                else
                    dist_f=zeros(size(c_tu_int_c,1),1);
                    for n_2=1:size(c_tu_int_c,1)
                        dist_f(n_2)=sqrt((data_tu(end,3)-(c_tu_int_c(n_2, 2)+x_c-50)).^2+(data_tu(end,2)-(c_tu_int_c(n_2,1)+y_c-50)).^2);
                    end
                    [~,I]=min(dist_f);
                end
                c_tu_int_c=c_tu_int_c(I,:);
                r_tu_int_c=r_tu_int_c(I);
            end
            clear dist_c final_idx I dist_f;
        end
        clear n_cells;
    end

    if exist('c_tu_int_c', 'var') && not(isempty(c_tu_int_c))
        c_tu_int_c = [c_tu_int_c(1,1)+y_c-50 c_tu_int_c(1,2)+x_c-50];
    end
    c_tu_c = c_tu_int_c;
    r_tu_c =r_tu_int_c;

    c_tu=c_tu_int;
    r_tu =r_tu_int;

    if exist('c_tu_c', 'var') && not(isempty(c_tu_c))
        data_tu=[data_tu;i c_tu_c(1,[1 , 2])];
    end

    % tracking immune cells
    if flag_imm
        [~ ,c_sp,r_sp]=my_CHT_cell(Ib2_sp,polarity,round(max(2,r_min_sp*fac-r_min_sp*fac*0.15)),round(max(6,r_min_sp*fac+r_min_sp*fac*0.15)),fac,0.95);
        % for each cancer cell identified, compute the distance
        % from this one and each immune cell.
        % if the distance < r_tu_int_c, delete the immune cell.
        if exist("c_tu", "var") && not(isempty(c_tu))
            %cancelled_sp=[];
            %cancelled_r_sp=[];
            for j_sp=1:size(c_tu)
                c_cell=c_tu(j_sp,:);
                dist_sp=zeros(size(c_sp,1),1);
                for i_sp=1:size(c_sp,1)
                    dist_sp(i_sp)=sqrt((c_sp(i_sp,2)-c_cell(1, 2)).^2+(c_sp(i_sp,1)-c_cell(1,1)).^2);
                end
                if exist("r_tu_int_c", "var") && not(isempty(r_tu_int_c))
                    sp_idx=dist_sp<r_tu_int_c;
                else
                    sp_idx=dist_sp<r_min_tu_int;
                end
                if any(sp_idx)
                    %cancelled_sp=[cancelled_sp; c_sp(sp_idx,:)];
                    %cancelled_r_sp=[cancelled_r_sp; r_sp(sp_idx)];
                    c_sp(sp_idx,:)=[];
                    r_sp(sp_idx,:)=[];                    
                end
                clear c_cell dist_sp sp_idx;
            end
        end

        if exist('c_sp', 'var') && not(isempty(c_sp))
            data_sp=[data_sp;repmat(i,length(r_sp),1) c_sp(:,[1 , 2])];
        end
    else
        c_sp=[];
        r_sp=[];
        %cancelled_sp=[];
        %cancelled_r_sp=[];
    end

    % show
    if mod(i,400)==0
        figure(1);
        imshow(Ib2_sp); hold on;
        viscircles(c_sp,r_sp,'EdgeColor','r');
        hold on; viscircles(c_tu_c,r_tu_c,'EdgeColor','b');
        % hold on; viscircles(cancelled_sp, cancelled_r_sp, 'EdgeColor', 'y');
        % hold on;
        % viscircles(c_tu,r_tu,'EdgeColor','g');
        drawnow;
        hold off;
        %pause()
    end


    clear Ib2_sp roi...
        c_tu_int_c r_tu_int_c c_tu_int r_tu_int ...
        c_tu_c r_tu_c c_tu r_tu...
        c_sp r_sp...
        dist_sp
end

%% sp
if flag_imm
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
%% tu
id_tu=data_tu(1,1);
if id_tu > 1
    id_d=data_tu(data_tu(:,1)==id_tu,:);
    for i = 1 : id_tu-1
        data_tu2((i-1)*size(id_d,1)+1:i*size(id_d,1),2:3)=data_tu(1:size(id_d,1),2:3);
        data_tu2((i-1)*size(id_d,1)+1:i*size(id_d,1),1)=i;
    end
    data_tu2(size(data_tu2,1)+1:size(data_tu2,1)+size(data_tu,1),:)=data_tu;
    clear data_tu
    data_tu=data_tu2;
    clear data_tu2 id_tu id_d;
end

id_tu_fin=data_tu(end,1);
id_pos=find(data_tu(:,1)==id_tu_fin, 1, 'last');
if id_tu_fin < Nframes
    id_d=data_tu(data_tu(:,1)==id_tu_fin,:);
    for i = 1:(Nframes-id_tu_fin)
        data_tu2((i-1)*size(id_d,1)+1:i*size(id_d,1),2:3)=data_tu(id_pos:(id_pos+size(id_d,1)-1),2:3);
        data_tu2((i-1)*size(id_d,1)+1:i*size(id_d,1),1)=i+id_tu_fin;
    end
    data_tu=[data_tu; data_tu2];
    clear data_tu2
    clear data_tu2 id_tu id_tu_fin;
end


end

% function to detect cells
function [CHT,centers,radii]=my_CHT_cell(B,polarity,rmin,rmax,fac,sens)
warning off;
% fac=1;
if fac>1
    I2=imresize(B,fac,'bilinear');
    %I2=my_morph_sharp(I2,round(rmin/2));
else
    I2=B;
end
centers =[];
% sens=0.96;%0.95;
[centers,radii,~]=imfindcircles(I2,[rmin rmax],'ObjectPolarity',polarity,'Sensitivity',sens);
while isempty(centers) && sens<0.99
    sens=min(0.99,sens+0.01);
    [centers,radii,~]=imfindcircles(I2,[rmin rmax],'ObjectPolarity',polarity,'Sensitivity',sens);
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


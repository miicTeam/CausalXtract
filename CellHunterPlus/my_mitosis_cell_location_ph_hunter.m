function [data_tu]=my_mitosis_cell_location_ph_hunter(vidFrames, indexes, r_min_tu_int, polarity, dist_tu)

close all
fac=1;
m=1;
Nframes =indexes(end);
data_tu=[];

r_min_tu_int=r_min_tu_int/m;
x_c=size(vidFrames,1)/2;
y_c=size(vidFrames,2)/2;
for i = indexes
    [i Nframes]
    Ib2_sp=im2double(vidFrames(:,:,i));
    Ib2_sp(Ib2_sp>1)=1;
    Ib2_sp(Ib2_sp<0)=0;
    % Ib2_sp=imadjust(Ib2_sp,stretchlim(Ib2_sp));
    roi=imcrop(Ib2_sp,[y_c-50 x_c-50 100 100]);
    % tracking tumor cells (two cells into the center of the frame)
    [~,c_tu_int_c,r_tu_int_c, metric]=my_CHT_cell(roi,polarity,round(max(6/m,r_min_tu_int*fac-r_min_tu_int*fac*0.10)),round(min(16/m,r_min_tu_int*fac+r_min_tu_int*fac*0.10)),fac,0.98);

    if exist("c_tu_int_c", "var") && ~isempty(c_tu_int_c)
        n_cells=size(c_tu_int_c,1);
        % if there is more than one cell (as it should be):
        % compute distances between the multiple cells and the center of
        % the roi;
        % if the distance per cell is less than dist_tu px, keep the cells;
        % if there are >2 cells, keep the ones with the highest metric;
        if n_cells>1
            dist_c=zeros(n_cells,1);
            %final_idx=zeros(n_cells,1);
            for num_c=1:n_cells
                dist_c(num_c)=sqrt((c_tu_int_c(num_c,2)-size(roi,1)/2).^2+((c_tu_int_c(num_c,1)-size(roi,2)/2).^2));
            end
            final_idx=dist_c<dist_tu;
            dist_c=dist_c(final_idx,:);
            c_tu_int_c=c_tu_int_c(final_idx,:);
            r_tu_int_c=r_tu_int_c(final_idx);
            metric=metric(final_idx);
            if size(dist_c,1)>2
                [~,I] = maxk(metric,2);
                c_tu_int_c=c_tu_int_c(I,:);
                r_tu_int_c=r_tu_int_c(I);
            end
            clear dist_c final_idx I metric;
        end
        clear n_cells;
    end

    if exist('c_tu_int_c', 'var') && not(isempty(c_tu_int_c))
        c_tu_int_c=[c_tu_int_c(:,1)+y_c-50 c_tu_int_c(:,2)+x_c-50];
        c_tu_c=[mean(c_tu_int_c(:,1),1) mean(c_tu_int_c(:,2),1)];
    end

    if exist('r_tu_int_c', "var") && not(isempty(r_tu_int_c))
        r_tu_c = mean(r_tu_int_c);
    end

    if exist('c_tu_c', "var") && not(isempty(c_tu_c))
        data_tu=[data_tu;i c_tu_c(1,[1 , 2])];
    end
    %
    if mod(i,400)==0
        figure(1);
        imshow(Ib2_sp); hold on;
        %
        hold on;
        if exist('c_tu_int_c', 'var') && not(isempty(r_tu_int_c))
            viscircles(c_tu_int_c, r_tu_int_c,'EdgeColor','b'); hold on;
            plot(c_tu_c(:,1), c_tu_c(:,2), '*', 'Color', [1 1 1]);
        end
        drawnow;
        title(i)
        hold off;
        pause(0);
    end
    clear Ib2_sp roi c_tu_int_c r_tu_int_c...
        c_tu_c r_tu_c;
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
end

function [CHT,centers,radii, metric]=my_CHT_cell(B,polarity,rmin,rmax,fac,sens)
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
[centers,radii,metric]=imfindcircles(I2,[rmin rmax],'ObjectPolarity',polarity,'Sensitivity',sens);
edge_th=graythresh(I2);
while all([size(centers,1)<2 sens<=0.99 edge_th>0])
    sens=min(0.99,sens+0.01);
    [centers,radii,metric]=imfindcircles(I2,[rmin rmax],'ObjectPolarity',polarity,'Sensitivity',sens, 'EdgeThreshold', edge_th);
    if all([size(centers,1)<2 sens==0.99])
        edge_th=edge_th-0.01;
        fin_edge_th=max(0.01, edge_th);
        [centers,radii,metric]=imfindcircles(I2,[rmin rmax],'ObjectPolarity',polarity,'Sensitivity',sens, 'EdgeThreshold', fin_edge_th);
    end
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




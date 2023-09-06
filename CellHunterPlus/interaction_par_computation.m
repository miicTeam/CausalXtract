function [count_cell_at_frame, count_cell_until_fr,sum_cum_int,max_cum_int, min_cum_int, min_dist,max_dist,mean_dist,IDInt] = interaction_par_computation(p_sp2,p_tu2,r_int)

d_tu0=r_int;
p_sp6=[p_sp2];
p_tu6=[p_tu2];
Nframe_tu= length([p_tu6.t]);
for i = 1 : length(p_sp6)
    L(i) = length(p_sp6(i).t);
end
Nframe_sp = max(L);
X_tu=nan(Nframe_tu,length(p_tu6));
Y_tu=nan(Nframe_tu,length(p_tu6));

if Nframe_sp < Nframe_tu
    X_sp=nan(Nframe_tu,length(p_sp6));
    Y_sp=nan(Nframe_tu,length(p_sp6));
else
    X_sp=nan(Nframe_sp,length(p_sp6));
    Y_sp=nan(Nframe_sp,length(p_sp6));
end
for j = 1 : length(p_tu6.t)
    X_tu(p_tu6.t(j),1)=p_tu6.x(j);
    Y_tu(p_tu6.t(j),1)=p_tu6.y(j);
end


for i = 1 : length(p_sp6)
    for j = 1 : length(p_sp6(i).t)
        X_sp(p_sp6(i).t(j),i)=p_sp6(i).x(j);
        Y_sp(p_sp6(i).t(j),i)=p_sp6(i).y(j);
    end
end

id_int0=zeros(length(p_sp6),Nframe_tu);
id_dist = id_int0;
for k = 1 : Nframe_tu
    k
    %     [i k]
    a=pdist2([X_tu(k) Y_tu(k)],[X_sp(k,:)' Y_sp(k,:)']);
    id=find(a <= d_tu0);
    if not(isempty(id))
        id_int0(id,k)=1;
        id_dist(id,k)=a(id);

    end
    count_cell_at_frame(k) = sum(id_int0(:,k));
    if count_cell_at_frame(k)>0
        IDInt(k) = 1;
    else
        IDInt(k) = 0;
    end
    [r,~]=find(id_int0(:,1:k)==1);
    if isempty(r)==1
        count_cell_until_fr(k)=0;
    else
        count_cell_until_fr(k) = sum(length(unique(r)));
    end
    sum_cum_int(k) = sum(sum(id_int0(:,1:k)));
    sum_cum = sum(id_int0(:,1:k),2);
    sum_cum(sum_cum==0)=[];
    if isempty(sum_cum)==1
        max_cum_int(k) = 0;
        min_cum_int(k) = 0;
    else
        max_cum_int(k) = max(sum_cum);
        min_cum_int(k) = min(sum_cum);
    end
    col_dist = id_dist(:,k);
    col_dist(col_dist==0)=[];
    if isempty(col_dist)==1
        min_dist(k) = 0;
        max_dist(k) = 0;
        mean_dist(k) = 0;
    else
        min_dist(k) = min(col_dist);
        max_dist(k) = max(col_dist);
        mean_dist(k) = mean(col_dist);
    end
    clear a id r sum_cum col_dist;
end
min_dist(find(IDInt==0))=nan;
max_dist(find(IDInt==0))=nan;
mean_dist(find(IDInt==0))=nan;
max_cum_int(find(IDInt==0)) = nan;
min_cum_int(find(IDInt==0)) = nan;
sum_cum_int(find(IDInt==0)) = nan;
count_cell_until_fr(find(IDInt==0)) = nan;
count_cell_at_frame(find(IDInt==0)) = nan;


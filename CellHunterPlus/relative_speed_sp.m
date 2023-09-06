function [mean_v,max_v,min_v] = relative_speed_sp(p_tu2,p_sp2,r_int)

d_tu0=r_int;

p_sp6=[p_sp2];
p_tu6=[p_tu2];
Nframe_tu= length([p_tu6.t]);

X_tu=nan(Nframe_tu,length(p_tu6));
Y_tu=nan(Nframe_tu,length(p_tu6));

%if Nframe_sp < Nframe_tu
X_sp=nan(Nframe_tu,length(p_sp6));
Y_sp=nan(Nframe_tu,length(p_sp6));
V_sp=nan(Nframe_tu,length(p_sp6));
% else
%     X_sp=nan(Nframe_sp,length(p_sp6));
%     Y_sp=nan(Nframe_sp,length(p_sp6));
%     V_sp=nan(Nframe_sp,length(p_sp6));
%end
for j = 1 : length(p_tu6.t)
    X_tu(p_tu6.t(j),1)=p_tu6.x(j);
    Y_tu(p_tu6.t(j),1)=p_tu6.y(j);
end


for i = 1 : length(p_sp6)

    j1=1;
    while j1<=Nframe_tu && j1<=length(p_sp6(i).t)
        X_sp(p_sp6(i).t(j1),i)=p_sp6(i).x(j1);
        Y_sp(p_sp6(i).t(j1),i)=p_sp6(i).y(j1);
        j1=j1+1;
    end
    j2=2;
    while j2<=Nframe_tu && j2<=length(p_sp6(i).t)
        rel_pos_x_j = p_sp6(i).x(j2)- X_tu(j2);
        rel_pos_x_jmeno1 = p_sp6(i).x(j2-1)- X_tu(j2-1);
        rel_pos_y_j = p_sp6(i).y(j2)- Y_tu(j2);
        rel_pos_y_jmeno1 = p_sp6(i).y(j2-1)- Y_tu(j2-1);
        V_sp(p_sp6(i).t(j2),i)=pdist2([rel_pos_x_j rel_pos_y_j],[rel_pos_x_jmeno1 rel_pos_y_jmeno1])/(p_sp6(i).t(j2)-p_sp6(i).t(j2-1));
        j2=j2+1;
    end
end
id_int0=zeros(length(p_sp6),Nframe_tu);
id_v = id_int0;

for k = 1 : Nframe_tu
    k
    %     [i k]
    a=pdist2([X_tu(k) Y_tu(k)],[X_sp(k,:)' Y_sp(k,:)']);
    id=find(a <= d_tu0);
    if not(isempty(id))
        id_v(id,k)=V_sp(k,id);
    end
    col_v = id_v(:,k);
    col_v(col_v==0)=[];
    col_v(isnan(col_v))=[];
    if isempty(col_v)==1
        min_v(k) = nan;
        max_v(k) = nan;
        mean_v(k) = nan;
    else
        min_v(k) = min(col_v);
        max_v(k) = max(col_v);
        mean_v(k) = nanmean(col_v);
    end
    clear a id col_v
end



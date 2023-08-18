function [p_sp2,q_sp2]=my_track_refining(p_sp,q_sp,r_min,lmin,flag)
%% This function eliminates null data at the beginning of trajectories
for i = 1 : length(p_sp)
    if p_sp(i).y(1)==0

        id=find(p_sp(i).y==0);
        p_sp(i).y(id)=[];
        p_sp(i).x(id)=[];
        p_sp(i).t(id)=[];

    end
end
p_sp2=p_sp;

if not(isempty(q_sp))
    for i = 1 : length(q_sp)
        if q_sp(i).y(1)==0

            id=find(q_sp(i).y==0);
            q_sp(i).y(id)=[];
            q_sp(i).x(id)=[];
            q_sp(i).t(id)=[];

        end
    end
end
q_sp2=q_sp;


% Delete stopping cells (not assigned)
clear p_sp2;
p_sp0=p_sp;
for l = 1 : length(p_sp)
        if strcmp(flag,'sp')
            p_sp=my_red_track(p_sp,2*r_min,lmin);
        else
            p_sp=my_refining_length(p_sp,lmin);
        end
end
p_sp2=p_sp;
p_sp=p_sp0;

% save([dir_name,nome_mat],'p_sp2','-append');
if not(isempty(q_sp))
    clear q_sp2;
    q_sp0=q_sp;
    for l = 1 : length(q_sp)
        if strcmp(flag,'sp')
            q_sp=my_red_track(q_sp,2*r_min,lmin);
        else
            q_sp=my_refining_length(q_sp,lmin);
        end
    end
    q_sp2=q_sp;
    q_sp=q_sp0;
else
    q_sp2=[];
end
end
function p_sp2=my_refining_length(p_sp,L)
for i = 1 : length(p_sp)
  if p_sp(i).y(1)==0
    
      id=find(p_sp(i).y==0);
      p_sp(i).y(id)=[];
      p_sp(i).x(id)=[];
      p_sp(i).t(id)=[];
  end
end
n=1;
for i = 1 : length(p_sp)
    if size(p_sp(i).x,2)<=L
    else
        p_sp2(n)=p_sp(i);
        n=n+1;
    end
end
if not(exist('p_sp2'))
p_sp2=[];
end
end

function p_sp2=my_red_track(p_sp,r_min,lmin)


id_kept=true(length(p_sp),1);
for i = 1 : length(p_sp)
    i;

if  std(sqrt((p_sp(i).y-p_sp(i).y(1)).^2+(p_sp(i).x-p_sp(i).x(1)).^2)) < round(r_min) || length(unique([p_sp(i).x' p_sp(i).y'],'rows')) < lmin
    id_kept(i)=false;
else
    id_kept(i)=true;
end
end
   p_sp2=p_sp(id_kept); 

end
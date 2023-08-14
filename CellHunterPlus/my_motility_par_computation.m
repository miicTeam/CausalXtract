function [IstSpeed, StraightIndex,TotDist,NetDist]= my_motility_par_computation(track,trackt,Ml)
dist_total = 0;
for i = 2:Ml
    vn1 = track(i,:);
    vn2 = track(i-1,:);
    v1 = track(1,:);
    DT = trackt(i)-trackt(i-1);
    IstSpeed(i) = pdist([vn2;vn1])/DT;
    DT1 = trackt(i)-trackt(1);
    NetDist(i)= pdist([vn1;v1]);
    dist_total = dist_total + IstSpeed(i)*DT;
    TotDist(i) = dist_total;
    StraightIndex(i) = NetDist(i) /TotDist(i);
end
IstSpeed(1) = nan;
StraightIndex(1) =  nan;
TotDist(1)=nan;
NetDist(1)=nan;

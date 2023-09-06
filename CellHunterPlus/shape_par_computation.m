function [Area, Perimeter,Circularity,Eccentricity,InstChangeShape]= shape_par_computation(roi_tu,Nframes,m)
vidFrames_tu = roi_tu;
i = 1;
while i<=Nframes
    disp(i)
    Ib2_tu=vidFrames_tu(:,:,i);
    mask = zeros(size(Ib2_tu));
    mask(floor(24/m):end-floor(24/m),floor(24/m):end-floor(24/m)) = 1;
    Ib3_tu = activecontour(Ib2_tu,mask,300);
    BW = logical( Ib3_tu );
    R=regionprops(BW,'Area','Perimeter','Circularity','Eccentricity','EquivDiameter');
    [max_a,idx_a]= max([R.Area]);
    Area(i) = max_a*m^2;
    per = [R.Perimeter];
    Perimeter(i) = per(idx_a)*m;
    cir =[R.Circularity];
    Circularity(i) = cir(idx_a);
    ecc = [R.Eccentricity];
    Eccentricity(i) = ecc(idx_a);
    eq_diam = [R.EquivDiameter];
    Equivdiameter(i) = eq_diam(idx_a);

    i=i+1;
    clear Ib2_tu BW mask R Ib3_tu max_a idx_a per cir ecc eq_diam
end

InstChangeShape0 = abs(diff(Equivdiameter));
InstChangeShape = [nan,InstChangeShape0]';
end


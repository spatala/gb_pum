function [C1, r2, TR1] = identify_gb_vv(n_pts, gb_inds, TR, C, r, lat_par)
i1 = zeros(n_pts,1);
i1(gb_inds) = 1;

%%% Adjust this for more Voronoi Vertices
rCut = lat_par;

%%%% Remove large "raddii" 
C1 = C(sum(i1(TR),2) > 0,:);
Idx = rangesearch(C,C1,rCut);
n_cc = size(C,1);
inds2 = zeros(100*n_cc,1);
i1 = 1;
for ct1 = 1:size(Idx,1)
    n1=size(Idx{ct1},2);
    i2 = i1 + n1 -1;
    inds2(i1:i2) = Idx{ct1};
    i1 = i2 + 1;
end
if (i1 < 100*n_cc)
    inds2(i1:end) = [];
end
inds2 = unique(inds2);

C1 = C(inds2,:); r1 = r(inds2); TR1 = TR(inds2,:);
r2 = r1 - lat_par/(2*sqrt(2));

end
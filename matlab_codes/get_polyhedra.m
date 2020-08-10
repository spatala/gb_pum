function [n_poly_atms, poly_inds, inpoly_inds] = get_polyhedra(C1,r2, TR1, pts, box_cell)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C1: The coordinates of Voronoi vertices
% r2: The radii of the Voronoi spheres (circum_sphere - atom)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cc_tet = delaunayTriangulation(C1);
cc_TR = cc_tet.ConnectivityList;

inds = [[1,2];[1,3];[1,4];[2,3];[2,4];[3,4]];

conn_nodes = zeros(size(cc_TR,1)*6,2);
ist = 1;
for ct1=1:6
    i1 = inds(ct1,1);
    i2 = inds(ct1,2);
    ti1 = cc_TR(:,i1);
    ti2 = cc_TR(:,i2);
    dpts = C1(ti1,:) - C1(ti2,:);
    d1 = sqrt(sum(dpts.^2,2)) - (r2(ti1) + r2(ti2));
    inpoly_inds = (d1 < 0);
    ti3 = ti1(inpoly_inds);
    ti4 = ti2(inpoly_inds);
    ni = size(ti3,1);
    istop = ist + ni-1;
    conn_nodes(ist:istop,:) = [ti3, ti4];
    ist = istop + 1;
end
conn_nodes(ist:end,:)=[];

unq_edges = unique(conn_nodes,'rows');
unq_edges = sort(unq_edges,2);
unq_edges = unique(unq_edges,'rows');

n_cc = size(C1,1);
A = sparse(n_cc,n_cc);
i1 = (unq_edges(:,1)-1)*(n_cc) + unq_edges(:,2);
A(i1) = 1;
i1 = (unq_edges(:,2)-1)*(n_cc) + unq_edges(:,1);
A(i1) = 1;

G1 = graph(A);
[bins,binsizes] = conncomp(G1);

n_polys = size(binsizes,2);
poly_inds = cell(n_polys,1);
n_poly_atms = zeros(n_polys,1);
cm_poly = zeros(n_polys,3);
for ct1=1:n_polys
    i1 = find(bins == ct1);
    i2 = TR1(i1,:);
    i3 = unique(i2(:));
    poly_inds{ct1} = i3;
    n_poly_atms(ct1) = size(i3,1);
    cm_poly(ct1,:) = sum(pts(i3,:))/size(i3,1);
end

inpoly_inds = inpoly_units(box_cell, cm_poly);


end
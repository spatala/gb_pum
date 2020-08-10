function [] = cc_coors_analysis(mat_name)
%% Load Data
s1 = load(mat_name);
% The grain-boundary atoms + atoms 
pts = s1.pts;
n_pts = size(pts,1);

gb_inds = s1.gb_inds+1;
uc_inds = s1.uc_inds+1;

tet = delaunayTriangulation(pts);
TR = tet.ConnectivityList;
[C,r] = circumcenter(tet);

lat_par = 4.05;

[C1, r2, TR1] = identify_gb_vv(n_pts, gb_inds, TR, C, r, lat_par);


%% get_polyhedra
box_cell = s1.box_cell;
[n_poly_atms, poly_inds, inpoly_inds] = ...
    get_polyhedra(C1,r2, TR1, pts, box_cell);

%% get_polyhedra_attributes
poly_attr = get_polyhedra_attributes(n_poly_atms, poly_inds, ...
    inpoly_inds, uc_inds);

%% find_sim_box_poly
poly_attr = find_sim_box_poly(poly_attr);


%% Plot polyhedral Units
plot_gb_poly(pts, gb_inds,poly_attr, mat_name)
end
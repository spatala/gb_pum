function poly_attr = get_polyhedra_attributes(n_poly_atms, poly_inds, ...
    inpoly_inds, uc_inds)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Get polyhedral units attributes
n_polys = size(inpoly_inds,1);
n_poly_atms = n_poly_atms(inpoly_inds);
poly_szs = unique(n_poly_atms);
poly_attr = cell(size(unique(n_poly_atms)));

ind_map = zeros(max(poly_szs),1);
ind_map(poly_szs) = 1:size(poly_szs,1);

for ct1 = 1:size(poly_szs,1)
    n1=poly_szs(ct1);
    n_poly = size(find(n_poly_atms == n1),1);
    
    poly_attr_n.n_atm = n1;
    poly_attr_n.atm_inds = zeros(n_poly,n1);
    poly_attr_n.uc_inds = zeros(n_poly,n1);
    
    poly_attr{ct1} = poly_attr_n;
end

n_ctr = zeros(max(poly_szs),1);
for ct1 = 1:n_polys
    n1 = n_poly_atms(ct1);
    atm_inds1 = poly_inds{inpoly_inds(ct1)};
    uc_inds1 = uc_inds(atm_inds1);
    poly_attr{ind_map(n1)}.atm_inds(n_ctr(n1)+1,:) = atm_inds1;
    poly_attr{ind_map(n1)}.uc_inds(n_ctr(n1)+1,:) = uc_inds1;
    n_ctr(n1) = n_ctr(n1) + 1;
end

end
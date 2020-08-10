function poly_attr = find_sim_box_poly(poly_attr)

for ct1=1:size(poly_attr,1)
    [~, ia] = unique(poly_attr{ct1}.uc_inds,'rows');
    atm_inds1 = poly_attr{ct1}.atm_inds(ia,:);
    uc_inds3 = poly_attr{ct1}.uc_inds(ia,:);
    poly_attr{ct1}.atm_inds = atm_inds1;
    poly_attr{ct1}.uc_inds = uc_inds3;
end

end
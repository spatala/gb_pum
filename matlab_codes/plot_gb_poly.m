function [] = plot_gb_poly(pts, gb_inds,poly_attr,mat_name)
%% Plot polyhedra units
figure('Position',[10,10,1200,1200]); hold on;

tpts = pts(gb_inds,:);
x1 = tpts(:,1); y1 = tpts(:,2); z1 = tpts(:,3);
scatter3(x1,y1,z1,1000,'.');

pol_color = {[226 212 255]/255, [13 188 68]/255, [107 139 240]/255,...
    [64,224,208]/255, [255,255,51]/255, [250 42 27]/255, [32 254 43]/255,...
    [58 96 33]/255, [98 105 85]/255, [96 53 89]/255, [175 53 205]/255, [0 255 255]/255,...
    [138,83,226]/255, [230 230 42]/255, [83,83,83]/255};

for ct1=1:size(poly_attr,1)
    n1 = poly_attr{ct1}.n_atm;
    if n1 > 6
        atm_inds1 = poly_attr{ct1}.atm_inds;
        if n1 < 18
            pol_col = pol_color{n1-3};
        else
            pol_col = [0.5, 0.5, 0.5];
        end
        n_pols = size(atm_inds1,1);
        for ct2 = 1:n_pols
            poly_i = atm_inds1(ct2,:);
            poly_pts = pts(poly_i,:);
            
            tet = delaunayTriangulation(poly_pts);
            [T,Xb] = freeBoundary(tet);
            TR = triangulation(T,Xb);
            trisurf(TR, poly_pts,'FaceColor', pol_col,...
                'FaceAlpha', 0.2, 'LineWidth',2,...
                'EdgeColor',pol_col,'EdgeAlpha',0.7);hold on;
        end
    end
end

zlim([-30,30]);

fig_name = [mat_name,'_gb_poly'];
export_fig(fig_name, '-a1', '-r200','-transparent','-png');

end
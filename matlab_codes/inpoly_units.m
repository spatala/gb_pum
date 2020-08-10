function ind1 = inpoly_units(box_cell, cm_poly)

xvec = box_cell(:,1);
yvec = box_cell(:,2);
zvec = box_cell(:,3);
orig = box_cell(:,4);
box_pts = zeros(8,3);
box_pts(1,:) = orig;
box_pts(2,:) = orig + xvec;
box_pts(3,:) = orig + xvec + yvec;
box_pts(4,:) = orig + yvec;
box_pts(5,:) = orig + zvec;
box_pts(6,:) = orig + zvec + xvec;
box_pts(7,:) = orig + zvec + xvec + yvec;
box_pts(8,:) = orig + zvec + yvec;

tet = delaunayTriangulation(box_pts);
[T,Xb] = freeBoundary(tet);


FV.vertices = Xb; 
FV.faces = T;
in = inpolyhedron(FV,cm_poly); 
ind1 = find(in);
end
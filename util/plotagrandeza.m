
function plotagrandeza(coord,conec, U)
patch('Vertices',coord,'Faces',conec,'FaceVertexCData',U,'FaceColor','interp'); 
shading interp
lighting phong
colormap jet
% view(3);
% axis equal
% axis vis3d
colorbar
% grid
end
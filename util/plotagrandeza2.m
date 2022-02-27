function plotagrandeza2(x,y,U)

%     figure
%     mesh(x,y,U)
    surf(x,y,U)
%     pcolor(x,y,U)
%     grid
    colormap jet
    xlabel('L_x')
    ylabel('L_y')
    axis('square')
    axis vis3d
    shading interp
    colorbar
end
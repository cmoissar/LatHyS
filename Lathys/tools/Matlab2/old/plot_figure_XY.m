hax = axes('Units','normalized','position',[.1  .1  .8  .6])
pcolor(field_plane_XY');
shading flat;
%set(gca,'DataAspectRatio',[1 1 1]);
 % Set color-limits
   c_min = -0.2;
   c_max = 0.2;
   c_limit = [c_min c_max];
 %  set(gca, 'CLim', c_limit);
   
   % scaler at side
   h_bar = colorbar('vert');
   scale_data = linspace( c_limit(1), c_limit(2), length( colormap ));
   set(hhax,_bar, 'CLim', c_limit);
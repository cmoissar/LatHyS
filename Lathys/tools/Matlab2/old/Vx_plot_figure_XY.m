%%%%%%%%%%%%%%%%%%%%%%%%%
% Vx_plot_figure
%------------------------
% This routine plot the Vx field
% in XY  plane
%
% R. Modolo
% UVSQ / LATMOS 
% Mars 2011
%%%%%%%%%%%%%%%%%%%%%%%%%
plane_XY = 'XY';
field_val_XY = 'Vx';
field_XY = double(Vx);
val_plane_XY = round(nc(3)/2.);
cell = nc(3)-1;
field_plane_XY  = zeros(nc(1),nc(2));
field_plane_XY(:,:) =  field_XY(:,:,val_plane_XY);

scrsz = get(0,'ScreenSize');

x0 = 20.;
y0 = 50.;
x1 = 1.*scrsz(3)/3;
y1 = 2.*scrsz(4)/3.;


 h_f1 = figure('Position',[x0 y0 x1 y1],'Name','Vx XY',...
     'Numbertitle','off');
 
 pos_x = 0.1;
 pos_y = 0.1;
 if (nc(1) > nc(2))
   size_x = 0.75;
   size_y = size_x*nc(2)/nc(1);
 else
     size_y = 0.75;
     size_x = size_y*nc(1)/nc(2);
 end
     
 hax = axes('Units','normalized','position',[pos_x pos_y size_x size_y]);
 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% color button
    x0 = 20;
    y0 = 550;
    x1 = 50;
    y1 = 20;
    uicontrol('Style', 'popup',...
           'String', 'jet|hsv|hot|cool|gray',...
           'Position', [x0 y0 x1 y1],...
           'Callback', @setmap);       % Popup function handle callback
                                       % Implemented as a subfunction
   uicontrol('Style','text',...
        'Position',[x0 y0+y1 x1 20],...
        'String','Color')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save button
uicontrol('Style', 'pushbutton', 'String', 'Save',...
        'Position', [20 20 50 20],...
        'Callback', @save_figure);        % Pushbutton string callback
                                   % that calls a MATLAB function
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % slider horizontal
     x0 = 100;
     y0 = 550;
     x1 = 250;
     y1 = 20;
 uicontrol('Style', 'slider',...
         'Min',1,'Max',cell,'SliderStep',[1./(cell-1) 0.8],'Value',val_plane_XY,...
         'Position', [x0 y0 x1 y1],...
         'Callback', {@plane_value,field_XY,field_plane_XY,plane_XY,field_val_XY,cell,val_plane_XY});   % Uses cell array function handle callback
                                     % Implemented as a subfunction with an argument
 uicontrol('Style','text',...
         'Position',[x0 y0+y1 x1 20],...
         'String','Plane value')
%plot_figure_XY;



pcolor(field_plane_XY');
 shading flat;
% set(gca,'DataAspectRatio',[1 1 1],'DataAspectRatioMode','Manual');
%axis square;
  % Set color-limits
    c_min = -0.2;
    c_max = 0.2;
    c_limit = [c_min c_max];
    set(gca, 'CLim', c_limit);
    
    % scaler at side
   % h_bar = colorbar('location','EastOutside');
   % scale_data = linspace( c_limit(1), c_limit(2), length( colormap ));
   % set(h_bar, 'CLim', c_limit); 

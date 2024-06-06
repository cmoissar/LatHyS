%%%%%%%%%%%%%%%%%%%%%%%%%
% Bz_plot_figure
%------------------------
% This routine plot the Bz field
% in YZ plane
%
% R. Modolo
% UVSQ / LATMOS 
% Mars 2011
%%%%%%%%%%%%%%%%%%%%%%%%%
plane_YZ = 'YZ';
field_val_YZ = 'Bz';
field_YZ = double(Bz);
val_plane_YZ = round(nc(1)/2.);
cell = nc(2)-1;
field_plane_YZ  = zeros(nc(2),nc(3));
field_plane_YZ(:,:) =  field_YZ(val_plane_YZ,:,:);

scrsz = get(0,'ScreenSize');

 h_f3 = figure('Position',[2.*scrsz(3)/3-20 50 1.*scrsz(3)/3 2.*scrsz(4)/3.],'Name','Bz YZ',...
     'Numbertitle','off');
 
 pos_x = 0.1;
 pos_y = 0.1;
 if (nc(1) > nc(3))
   size_x = 0.7;
   size_y = size_x*nc(3)/nc(2);
 else
     size_y = 0.7;
     size_x = size_y*nc(2)/nc(3);
 end
     
 hax = axes('Units','normalized','position',[pos_x pos_y size_x size_y]);
 pcolor(field_plane_YZ');
 shading flat;
 set(gca,'DataAspectRatio',[1 1 1]);
  % Set color-limits
    c_min = -0.2;
    c_max = 0.2;
    c_limit = [c_min c_max];
    set(gca, 'CLim', c_limit);
    
    % scaler at side
   % h_bar = colorbar('location','EastOutside');
   % scale_data = linspace( c_limit(1), c_limit(2), length( colormap ));
   % set(h_bar, 'CLim', c_limit); 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% color button
    x0 = 20;
    y0 = 2.*scrsz(4)/3.-120;
    x1 = 50;
    y1 = 100;
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
     x0 = 300;
     y0 = 2.*scrsz(4)/3.-45.;
     x1 = 250;
     y1 = 20;
 uicontrol('Style', 'slider',...
         'Min',1,'Max',cell,'SliderStep',[1./(cell-1) 0.8],'Value',val_plane_YZ,...
         'Position', [x0 y0 x1 y1],...
         'Callback', {@plane_value,field_YZ,field_plane_YZ,plane_YZ,field_val_YZ,cell,val_plane_YZ});   % Uses cell array function handle callback
                                     % Implemented as a subfunction with an argument
 uicontrol('Style','text',...
         'Position',[x0 y0+y1 x1 20],...
         'String','Plane value')


%%%%%%%%%%%%%%%%%%%%%%%%%
% R. Modolo
% UVSQ / LATMOS 
% Mars 2011
%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_3_planes(field,field_val,posy,centr,radius)
h=findall(gcf,'-not','Style','Text','-and','-not','Style','Radio','-and','-not','Type','Figure','-and','-not','Type','uipanel');
delete(h);
plot_figure(field,'XY',field_val,1,posy,centr,radius);
plot_figure(field,'YZ',field_val,2,posy,centr,radius);
plot_figure(field,'XZ',field_val,3,posy,centr,radius);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% color bar and button
   scrsz = get(0,'ScreenSize');
   pos_x = 0.95;
   pos_y = 0.05;
   size_x=0.5;
   size_y=0.75*0.3*scrsz(3)/scrsz(4)*1.5

   hax = axes('Units','normalized','visible','off',...
       'position',[pos_x pos_y size_x size_y]);

    y0 = 2.*scrsz(4)/3.;%-85.;
    y1 = 20;
    x0 = pos_x*scrsz(3);
    x1 = 50;
   uicontrol('Style', 'popup',...
           'String', 'jet|hsv|hot|cool|gray',...
           'Position', [x0 y0 x1 y1],...
           'Callback', @setmap);       % Popup function handle callback
                                       % Implemented as a subfunction
   uicontrol('Style','text',...
        'Position',[x0 y0+y1 x1 20],...
        'String','Color')

   c_limit=[min(min(min(field))) max(max(max(field)))]
   caxis(hax,c_limit); 
   h_bar = colorbar('location','Westoutside');
   scale_data = linspace( c_limit(1), c_limit(2), length( colormap ));
   set(h_bar, 'CLim', c_limit); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end


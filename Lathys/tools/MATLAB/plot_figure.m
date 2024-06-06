%%%%%%%%%%%%%%%%%%%%%%%%%
% R. Modolo
% UVSQ / LATMOS 
% Mars 2011
%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  !!! All graphics modification is this function should
%      also appear in the plane_value function in order
%      to persist when the plane is changed
%
%%%%%%%%%%%%%%%%%%%%%%%%%
function c=plot_figure(field,plane,field_val,posx,posy,centr,radius)
ncells=size(field);
nc=ncells*0.;
 if plane=='XY'
        nc=ncells;
        cell = nc(3)-1;
        val_plane = round(nc(3)/2.);
        field_plane(:,:) = field(:,:,val_plane);
        centrx=centr(1);centry=centr(2);centrz=centr(3);
    elseif plane == 'XZ'        
        nc(1)=ncells(1);nc(2)=ncells(3);nc(3)=ncells(2);
        cell = nc(3)-1;
        val_plane = round(nc(3)/2.);
        field_plane(:,:) = field(:,val_plane,:);
        centrx=centr(1);centry=centr(3);centrz=centr(2);
    elseif plane=='YZ'
        nc(1)=ncells(2);nc(2)=ncells(3);nc(3)=ncells(1);
        cell = nc(3)-1;
        val_plane = round(nc(3)/2.);
        field_plane(:,:) = field(val_plane,:,:);
        centrx=centr(2);centry=centr(3);centrz=centr(1);
    end
  

scrsz = get(0,'ScreenSize');

 pos_x = 0.03+(1.*posx-1.)*0.31;
 pos_y = 0.05;
 if (nc(1) > nc(2))
   size_x = 0.75;
   size_y = size_x*nc(2)/nc(1);
 else
     size_y = 0.75;
     size_x = size_y*nc(1)/nc(2);
 end

 if (size_x<0.75*size_y) 
        size_x=size_x*1.3;
 end
 
size_x=size_x*0.37;
size_y=size_y*0.37*scrsz(3)/scrsz(4)*1.2;
 hax = axes('Units','normalized','position',[pos_x pos_y size_x size_y]);

 
 y0 = 2.*scrsz(4)/3.-85.;
    y1 = 20;
    x0 = 20+(1.*posx-1.)*0.31*scrsz(3);
    x1 = 50;
   
% save button
    x0 = (1.*posx-1.)*0.31*scrsz(3)+scrsz(3)/6+50.;
    x1 = 50;
   uicontrol('Style', 'pushbutton', 'String', 'Save',...
        'Position', [x0 y0 x1 y1],...
        'Callback', @save_figure);        % Pushbutton string callback
                                   % that calls a MATLAB function
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % slider horizontal
    x0 = (1.*posx-1.)*0.31*scrsz(3)+20.;
    x1 = scrsz(3)/6;

    uicontrol('Style','text',...
         'Position',[x0 y0+y1 x1-30 20],...
         'String',['Plane ',plane,',   value:']);
    ctext=uicontrol('Style','text',...
         'Position',[x0+x1-30 y0+y1 30 20],...
         'String',int2str(val_plane));

     pcolor(hax,field_plane');
    shading flat;
    c_limit=[min(min(min(field))) max(max(max(field)))];
    %c_limit=[min(min(field_plane)) max(max(field_plane))]
    caxis(hax,c_limit);

    
    x_p=0.0:1.0:360.0;
    dz=abs(val_plane-centrz)/radius;
    if (dz > 1)
        rd=0.;
    else
        rd=radius*sin(acos(dz));
    end
    y_p=sin(x_p*0.0174533)*rd+centry+1.5;
    x_p=cos(x_p*0.0174533)*rd+centrx+1.5;
    hold all;
    plot(hax,x_p,y_p,'Color','white','LineWidth',2);
    
    uicontrol('Style', 'slider',...
         'Min',1,'Max',cell,'SliderStep',[1./(cell-1) 0.8],'Value',val_plane,...
         'Position', [x0 y0 x1 y1],...
         'Callback', {@plane_value,field,field_plane,plane,field_val,cell,val_plane,hax,ctext,centr,radius}); 
        % Uses cell array function handle callback
        % Implemented as a subfunction with an argument

end
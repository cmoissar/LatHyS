function plane_value(hObj,event,field,field_plane,plane,field_val,cell,val_plane,h,htex,centr,radius) %#ok<INUSL>
    % Called to set zlim of surface in figure axes
    % when user moves the slider control
    axes(h)
    val_plane = round(get(hObj,'Value'));
    cell3 = size(field);
    if plane=='XY'
        cell = cell3(3)-1;
        field_plane(:,:) = field(:,:,val_plane);
        centrx=centr(1);centry=centr(2);centrz=centr(3);
    elseif plane == 'XZ'
        cell = cell3(2)-1;
        field_plane(:,:) = field(:,val_plane,:);
        centrx=centr(1);centry=centr(3);centrz=centr(2);
    elseif plane=='YZ'
        cell = cell3(1)-1;
        field_plane(:,:) = field(val_plane,:,:);
        centrx=centr(2);centry=centr(3);centrz=centr(1);
    end
    delete(findall(gca,'Type','surface'));
    delete(findall(gca,'Type','line'));
    pcolor(h,field_plane');
    shading flat;
    c_limit=[min(min(min(field))) max(max(max(field)))];
    %c_limit=[min(min(field_plane)) max(max(field_plane))]
    caxis(h,c_limit);
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
    plot(h,x_p,y_p,'Color','white','LineWidth',2);
    
    
    scrsz = get(0,'ScreenSize');
    x0 = scrsz(3)/9;
    x1 = scrsz(3)/6;
    y0 = 2.*scrsz(4)/3.-45.;
    y1 = 20;

set(htex,'String',int2str(val_plane));
end
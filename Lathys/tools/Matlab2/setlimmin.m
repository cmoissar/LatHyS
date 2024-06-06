function setlimmin(hObj,event) %#ok<INUSD>
        % Called when user activates popup menu 
        val = get(hObj,'String');
        h = findall(gcf,'Type','axes');
        c = caxis(h(2));
        c(1) = eval(val);
        caxis(h(1),c); 
        caxis(h(2),c); 
        caxis(h(3),c); 
        caxis(h(4),c); 
        caxis(h(5),c); 

end
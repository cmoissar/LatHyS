function setfield(hObj,event,list_obj)
what_field=get(hObj,'Value')
lst_cmp_spe=get(list_obj{35},'String')
if ((what_field == 2) | (what_field == 6))
        list_compo={'x' 'y' 'z' '//' '_|_' 'total'};
end
if ((what_field == 5))
        list_compo={'x' 'y' 'z' '//' '_|_' 'Pedersen' 'Hall' 'total'};
end
if ((what_field == 3))
        list_compo={'x' 'y' 'z' 'total'};
end
if ((what_field == 1) | (what_field == 4))
        list_compo={''};
end
if ((what_field == 7))
        list_compo=lst_cmp_spe;%{'electrons'};
end
if ((what_field == 8))
        list_compo=lst_cmp_spe;%{'electrons'};
end
if ((what_field == 9))
        list_compo={''};
end
obj=list_obj{26}; %u22
set(obj,'String',list_compo);
set(obj,'Value',1);
plot_figs(0,0,list_obj);
end
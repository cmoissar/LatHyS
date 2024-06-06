%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cw_read_MPI_netCDF.m
%--------------------------
% This program reads the cw files
% of the simulation written in
% netdcf format.
% It plots  : density, speeed,
% magnetic field, electric field
%--------------------------
% R. Modolo, M. Mancini
% LATMOS / UVSQ
% Quartier des Garennes
% 11 bd d'Alembert
% 78280 Guyancourt
% ronan.modolo@latmos.ipsl.fr


clear all;

typefile = 'Magw_';
runname = '06_06_11_';
dirname = '/home/seb/Bureau/Hybrid/CDF/1/';
src_name=[dirname typefile runname '*.nc'];
list=(dir(src_name));
src_name=[typefile runname 't'];

a=regexprep(list(1).name, src_name, '');
a=regexprep(a, '.nc', '');
list_time=a;
diagtime=strcat('t',a)
for i=2:size(list)
a=regexprep(list(i).name, src_name, '');
a=regexprep(a, '.nc', '');
list_time=strcat(list_time,'|');
list_time=strcat(list_time,a);
end

ncfile = [dirname typefile runname diagtime '.nc'];

ncid = netcdf.open(ncfile,'NC_NOWRITE');
% Get information about the contents of the file.
[numdims, numvars, numglobalatts, unlimdimID] = netcdf.inq(ncid);
nom_variable ='' ;

for i=0:numvars-1

[varname, xtype, dimids, numatts] = netcdf.inqVar(ncid,i);

nom_variable = char(nom_variable, varname)

end

planetname = transpose(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'planetname')))
centr      = transpose(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'s_centr')))
radius     = transpose(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'r_planet')))*1.
gs=transpose(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'gstep')))
radius=radius/gs(1);
centr=centr./gs;
Bx         = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'Bfield_x'));
nc = size(Bx);
Bx=0;
netcdf.close(ncid);

scrsz = get(0,'ScreenSize');
h_f1 = figure('Position',scrsz,'Name',' ',...
    'Numbertitle','off');

x0=scrsz(1)+10.;
y0=scrsz(2)+2.*scrsz(4)/3.-40.;
x1=scrsz(3);
y1=scrsz(4);
disp('y0');
disp(y0);

% We define 3 window plots for 3 slices in XY / XZ / YZ
% 1st Window
% ========== XY =============
 %Create the button group.
q0 = uicontrol('Style','frame',...
    'pos',[x0 y0 560 190]);
q1 = uicontrol('Style','frame',...
    'pos',[x0+570. y0 250 190]);

 h = uibuttongroup('visible','off','Position',[0 0 1 1]);

t0 = uicontrol('Style','Text','String','Planet:',...
    'pos',[x0 y0+191 75 20],'parent',h,'HandleVisibility','off');
t1 = uicontrol('Style','Text','String',planetname,...
    'pos',[x0+80 y0+191 75 20],'parent',h,'HandleVisibility','off');
% Create three radio buttons in the button group.
u0 = uicontrol('Style','Radio','String','Ex',...
    'pos',[x0 y0+150 75 30],'parent',h,'HandleVisibility','off');
u1 = uicontrol('Style','Radio','String','Ey',...
    'pos',[x0 y0+105 75 30],'parent',h,'HandleVisibility','off');
u2 = uicontrol('Style','Radio','String','Ez',...
    'pos',[x0 y0+60 75 30],'parent',h,'HandleVisibility','off');
u3 = uicontrol('Style','Radio','String','Bx',...
    'pos',[x0+80. y0+150 75 30],'parent',h,'HandleVisibility','off');
u4 = uicontrol('Style','Radio','String','By',...
    'pos',[x0+80. y0+105 75 30],'parent',h,'HandleVisibility','off');
u5 = uicontrol('Style','Radio','String','Bz',...
    'pos',[x0+80. y0+60 75 30],'parent',h,'HandleVisibility','off');
u6 = uicontrol('Style','Radio','String','Btot',...
    'pos',[x0+80. y0+15 75 30],'parent',h,'HandleVisibility','off');
u7 = uicontrol('Style','Radio','String','Vx',...
    'pos',[x0+160. y0+150 75 30],'parent',h,'HandleVisibility','off');
u8 = uicontrol('Style','Radio','String','Vy',...
    'pos',[x0+160. y0+105 75 30],'parent',h,'HandleVisibility','off');
u9 = uicontrol('Style','Radio','String','Vz',...
    'pos',[x0+160. y0+60 75 30],'parent',h,'HandleVisibility','off');
u10 = uicontrol('Style','Radio','String','Vtot',...
    'pos',[x0+240. y0+60 75 30],'parent',h,'HandleVisibility','off');
u14 = uicontrol('Style','Radio','String','V//',...
    'pos',[x0+240. y0+150 75 30],'parent',h,'HandleVisibility','off');
u15 = uicontrol('Style','Radio','String','Vperp',...
    'pos',[x0+240. y0+105 75 30],'parent',h,'HandleVisibility','off');
u16 = uicontrol('Style','Radio','String','Jx',...
    'pos',[x0+320. y0+150 75 30],'parent',h,'HandleVisibility','off');
u17 = uicontrol('Style','Radio','String','Jy',...
    'pos',[x0+320. y0+105 75 30],'parent',h,'HandleVisibility','off');
u18 = uicontrol('Style','Radio','String','Jz',...
    'pos',[x0+320. y0+60 75 30],'parent',h,'HandleVisibility','off');
u19 = uicontrol('Style','Radio','String','J//',...
    'pos',[x0+400. y0+150 75 30],'parent',h,'HandleVisibility','off');
u20 = uicontrol('Style','Radio','String','Jperp',...
    'pos',[x0+400. y0+105 75 30],'parent',h,'HandleVisibility','off');
u21 = uicontrol('Style','Radio','String','Jtot',...
    'pos',[x0+400. y0+60 75 30],'parent',h,'HandleVisibility','off');
u11 = uicontrol('Style','Radio','String','Dn',...
    'pos',[x0+480. y0+150 75 30],'parent',h,'HandleVisibility','off');
u12 = uicontrol('Style','Radio','String','Log(Dn)',...
    'pos',[x0+480. y0+105 75 30],'parent',h,'HandleVisibility','off');
u13 = uicontrol('Style','Radio','String','Sat. Dn',...
    'pos',[x0+480. y0+60 75 30],'parent',h,'HandleVisibility','off');
%u22 = uicontrol('Style','Radio','String','Prod_O',...
%    'pos',[x0+560. y0+150 75 30],'parent',h,'HandleVisibility','off');
%u23 = uicontrol('Style','Radio','String','Prod_H',...
%    'pos',[x0+560. y0+105 75 30],'parent',h,'HandleVisibility','off');
%u24 = uicontrol('Style','Radio','String','Prod_CO2',...
%    'pos',[x0+560. y0+60 75 30],'parent',h,'HandleVisibility','off');
%u2 = uicontrol('Style','Radio','String','Option 3',...
%    'pos',[10 150 100 30],'parent',h,'HandleVisibility','off');
% Initialize some button group properties. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% time button
x0=scrsz(1)+10.;
y0=scrsz(2)+2.*scrsz(4)/3.-40.;
src_name=[dirname 'Magw_' runname '*.nc'];
list=dir(src_name);
src_name=['Magw_' runname 't'];

b=regexprep(diagtime, 't', '');
list_time=b;
for i=1:size(list)
a=regexprep(list(i).name, src_name, '');
a=regexprep(a, '.nc', '');
    if(~strcmp(a,b))
        list_time=strcat(list_time,'|');
        list_time=strcat(list_time,a);
    end
end
 h2= uicontrol('Style', 'popup',...
           'String', list_time,...
           'Position', [x0+400 y0+185 75 30],...
           'Callback', {@settime,h,10,0,runname,diagtime,dirname,centr,radius})      % Popup function handle callback
                                       % Implemented as a subfunction
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Atmospheric density
typefile = 'Atmw_';
ncfile = [dirname typefile runname diagtime '.nc'];
ncid = netcdf.open(ncfile,'NC_NOWRITE');
% Get information about the contents of the file.
[numdims, numvars, numglobalatts, unlimdimID] = netcdf.inq(ncid);
nom_variable ='' ;
j=0;
for i=0:numvars-1
[varname, xtype, dimids, numatts] = netcdf.inqVar(ncid,i)
nom_variable = char(nom_variable, varname)
end
list_D=strmatch('Den_',nom_variable)
for i=1:size(list_D)
        k=fix((i-1)/3);
        l=mod((i-1),3);
        uicontrol('Style','Checkbox','String',regexprep(nom_variable(list_D(i),:), 'Den_', ''),...
        'pos',[x0+600.+k*90. y0+150-l*45. 89 30],'HandleVisibility','off','Callback',...
        {@dummy,h2,h,10,0,runname,diagtime,dirname,centr,radius});
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(h,'SelectionChangeFcn',{@selcbk,nc,h_f1,runname,diagtime,dirname,centr,radius,h2});
set(h,'SelectedObject',[]);  % No selection
set(h,'Visible','on');


function selcbk(source,eventdata,nc,h,runname,diagtime,dirname,centr,radius,h2)

function x=dneu
typefile ='Atmw_';
ncfile = [dirname typefile runname diagtime '.nc'];
ncid = netcdf.open(ncfile,'NC_NOWRITE');
x=dna;
test=0;
q=findall(gcf,'Style','Checkbox')
    for i=1:size(q)
        if (get(q(i),'Value') == 1)
            if (test == 0)
                x=x.*0;
                test =  1;
            end
            q(i)
                nom=strcat('Den_',get(q(i),'String'));
                x=x+netcdf.getVar(ncid,netcdf.inqVarID(ncid,nom));
        end
    end
netcdf.close(ncid);        
end

function x=dna
typefile ='Denw_';
ncfile = [dirname typefile runname diagtime '.nc'];
ncid = netcdf.open(ncfile,'NC_NOWRITE');
nrm   = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'phys_density'))
x     = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'Dn_tot'))*1E-6;%*2.5*1E-6;
x=x*nrm;
netcdf.close(ncid)
end

function x=Btot
x     = sqrt(Bx.^2+By.^2+Bz.^2);
end

function x=Bx
typefile ='Magw_';
ncfile = [dirname typefile runname diagtime '.nc'];
ncid = netcdf.open(ncfile,'NC_NOWRITE');
x     = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'Bfield_x'));
netcdf.close(ncid)
end

function x=By
typefile ='Magw_';
ncfile = [dirname typefile runname diagtime '.nc'];
ncid = netcdf.open(ncfile,'NC_NOWRITE');
x     = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'Bfield_y'));
netcdf.close(ncid)
end

function x=Bz
typefile ='Magw_';
ncfile = [dirname typefile runname diagtime '.nc'];
ncid = netcdf.open(ncfile,'NC_NOWRITE');
x     = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'Bfield_z'));
netcdf.close(ncid)
end

function x=Vtot
x     = sqrt(Vx.^2+Vy.^2+Vz.^2);
end

function x=Vx
typefile ='Velw_';
ncfile = [dirname typefile runname diagtime '.nc'];
ncid = netcdf.open(ncfile,'NC_NOWRITE');
x     = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'Vbulk_x'));
%x=x./dna;
netcdf.close(ncid)
end

function x=Vy
typefile ='Velw_';
ncfile = [dirname typefile runname diagtime '.nc'];
ncid = netcdf.open(ncfile,'NC_NOWRITE');
x     = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'Vbulk_y'));
%x=x./dna;
netcdf.close(ncid)
end

function x=Vz
typefile ='Velw_';
ncfile = [dirname typefile runname diagtime '.nc'];
ncid = netcdf.open(ncfile,'NC_NOWRITE');
x     = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'Vbulk_z'));
%x=x./dna;
netcdf.close(ncid)
end

function x=Ex
typefile ='Elew_';
ncfile = [dirname typefile runname diagtime '.nc'];
ncid = netcdf.open(ncfile,'NC_NOWRITE');
x     = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'Efield_x'));
netcdf.close(ncid)
end

function x=Ey
typefile ='Elew_';
ncfile = [dirname typefile runname diagtime '.nc'];
ncid = netcdf.open(ncfile,'NC_NOWRITE');
x     = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'Efield_y'));
netcdf.close(ncid)
end

function x=Ez
typefile ='Elew_';
ncfile = [dirname typefile runname diagtime '.nc'];
ncid = netcdf.open(ncfile,'NC_NOWRITE');
x     = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'Efield_z'));
netcdf.close(ncid)
end

function x=Jx
typefile ='Magw_';
ncfile = [dirname typefile runname diagtime '.nc'];
ncid = netcdf.open(ncfile,'NC_NOWRITE');
B_x     = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'Bfield_x'));
B_y     = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'Bfield_y'));
B_z     = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'Bfield_z'));
[J_x,J_y,J_z]=curl(B_x,B_y,B_z);
x     = J_x;
netcdf.close(ncid)
end

function x=Jy
typefile ='Magw_';
ncfile = [dirname typefile runname diagtime '.nc'];
ncid = netcdf.open(ncfile,'NC_NOWRITE');
B_x     = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'Bfield_x'));
B_y     = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'Bfield_y'));
B_z     = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'Bfield_z'));
[J_x,J_y,J_z]=curl(B_x,B_y,B_z);
x     = J_y;
netcdf.close(ncid)
end

function x=Jz
typefile ='Magw_';
ncfile = [dirname typefile runname diagtime '.nc'];
ncid = netcdf.open(ncfile,'NC_NOWRITE');
B_x     = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'Bfield_x'));
B_y     = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'Bfield_y'));
B_z     = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'Bfield_z'));
[J_x,J_y,J_z]=curl(B_x,B_y,B_z);
x     = J_z;
netcdf.close(ncid)
end

function x=Jtot
typefile ='Magw_';
ncfile = [dirname typefile runname diagtime '.nc'];
ncid = netcdf.open(ncfile,'NC_NOWRITE');
B_x     = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'Bfield_x'));
B_y     = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'Bfield_y'));
B_z     = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'Bfield_z'));
[J_x,J_y,J_z]=curl(B_x,B_y,B_z);
x = sqrt(J_x.^2+J_y.^2+J_z.^2);
netcdf.close(ncid)
end

function x=Jpar
typefile ='Magw_';
ncfile = [dirname typefile runname diagtime '.nc'];
ncid = netcdf.open(ncfile,'NC_NOWRITE');
B_x     = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'Bfield_x'));
B_y     = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'Bfield_y'));
B_z     = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'Bfield_z'));
[J_x,J_y,J_z]=curl(B_x,B_y,B_z);
x = (J_x.*B_x+J_y.*B_y+J_z.*B_z)./sqrt(B_x.^2+B_y.^2+B_z.^2);
netcdf.close(ncid)
end

figure(h);

Obj=get(h2,'String');
val = get(h2,'Value');
val=Obj(val,:);
diagtime=strcat('t',val); 

value = get(eventdata.NewValue,'String');
switch value
    case 'Ex'
        plot_3_planes(double(Ex),value,1,runname,diagtime,dirname,centr,radius);
    case 'Ey'
        plot_3_planes(double(Ey),value,1,runname,diagtime,dirname,centr,radius);
    case 'Ez'
        plot_3_planes(double(Ez),value,1,runname,diagtime,dirname,centr,radius);
    case 'Bx'
        plot_3_planes(double(Bx),value,1,runname,diagtime,dirname,centr,radius);
    case 'By'
        plot_3_planes(double(By),value,1,runname,diagtime,dirname,centr,radius);
    case 'Bz'
        plot_3_planes(double(Bz),value,1,runname,diagtime,dirname,centr,radius);
    case 'Btot'
        plot_3_planes(double(Btot),value,1,runname,diagtime,dirname,centr,radius);
    case 'Vx'
        plot_3_planes(double(Vx),value,1,runname,diagtime,dirname,centr,radius);
    case 'Vy'
        plot_3_planes(double(Vy),value,1,runname,diagtime,dirname,centr,radius);
    case 'Vz'
        plot_3_planes(double(Vz),value,1,runname,diagtime,dirname,centr,radius);
    case 'Vtot' 
        plot_3_planes(double(Vtot),value,1,runname,diagtime,dirname,centr,radius);
    case 'V//' 
        plot_3_planes(double((Vx.*Bx+Vy.*By+Vz.*Bz)./Btot),value,1,runname,diagtime,dirname,centr,radius);
    case 'Vperp' 
        plot_3_planes(double(sqrt(Vtot.^2-((Vx.*Bx+Vy.*By+Vz.*Bz)./Btot).^2)),value,1,runname,diagtime,dirname,centr,radius);
    case 'Dn'
        plot_3_planes(double(dneu),value,1,runname,diagtime,dirname,centr,radius);
    case 'Log(Dn)'
        d=dneu;
        mn=min(min(min(d(2:nc(1)-2,2:nc(2)-2,2:nc(3)-2))));
        plot_3_planes(double(log10(max(dneu,max(mn,1e-3)))),value,1,runname,diagtime,dirname,centr,radius);
    case 'Sat. Dn'
        plot_3_planes(double(min(dneu,2.*median(median(median(dneu))))),value,1,runname,diagtime,dirname,centr,radius);
    case 'Jx'
        plot_3_planes(double(Jx),value,1,runname,diagtime,dirname,centr,radius);
    case 'Jy'
        plot_3_planes(double(Jy),value,1,runname,diagtime,dirname,centr,radius);
    case 'Jz'
        plot_3_planes(double(Jz),value,1,runname,diagtime,dirname,centr,radius);
    case 'J//' 
        plot_3_planes(double(Jpar),value,1,runname,diagtime,dirname,centr,radius);
    case 'Jperp' 
        plot_3_planes(double(sqrt(Jtot.^2-Jpar.^2)),value,1,runname,diagtime,dirname,centr,radius);
    case 'Jtot' 
        plot_3_planes(double(Jtot),value,1,runname,diagtime,dirname,centr,radius);
    otherwise
        disp('Unknown Quantity');

end
end
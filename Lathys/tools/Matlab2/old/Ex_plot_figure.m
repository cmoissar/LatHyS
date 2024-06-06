%%%%%%%%%%%%%%%%%%%%%%%%%
% Ex_plot_figure
%------------------------
% This routine plot the Ex field
% in XY , XZ and YZ plane
%
% R. Modolo
% UVSQ / LATMOS 
% Mars 2011
%%%%%%%%%%%%%%%%%%%%%%%%%
%Ex_plot_figure_XY;
%Ex_plot_figure_XZ;
%Ex_plot_figure_YZ;

plot_figure(double(Ex),'XY','Ex',1,1);
plot_figure(double(Ex),'XZ','Ex',2,1);
plot_figure(double(Ex),'YZ','Ex',3,1);


function plot3_inertial(x,y,z,units,color)
%PLOT3_INERTIAL plots in inertial frame
%   Uses plot3, but automates formatting for 3d orbits.

xlabel_raw = '$\hat{\bf n_{1}}$ (' + units + ')';
xlabel(xlabel_raw,'interpreter','latex')
ylabel_raw = '$\hat{\bf n_{2}}$ (' + units + ')';
ylabel(ylabel_raw,'interpreter','latex')
zlabel_raw = '$\hat{\bf n_{3}}$ (' + units + ')';
zlabel(zlabel_raw,'interpreter','latex')
plot3(x,y,z,"Color",color);
end


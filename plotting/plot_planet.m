function plot_planet(radius,center,color)
%PLOT_PLANET Plots a planet in 3d
%   This function makes a "planet", or a sphere, in 3D by using the
%   bodies' mean equitorial radius and cartesian center point. inputs:
%   radius: planet's mean equitorial radius
%   center: 3 element vector specifying cartesian center point
%   color: color of planet, can be MATLAB color name or RGB triplet

[X,Y,Z] = sphere(30);
X = X * radius + center(1);
Y = Y * radius + center(2);
Z = Z * radius + center(3);
surf(X,Y,Z,'FaceColor',color,'EdgeColor','none')
axis equal
end


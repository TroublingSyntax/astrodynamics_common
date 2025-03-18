function plot_polyhedrons(F,V,Color,EdgeColor,EdgeThickness)
%PLOT_POLYHEDRONS Plots a general list of polyhedrons
%   This function will take on a list of vertices and faces and uses the
%   built-in patch() function in MATLAB to plot a 3D representation of the
%   model
%   INPUTS
%       F: the matrix of faces that connects the vertex indices listed in
%          each row
%       V: the matrix of vertices with each row being an (X,Y,Z) pair
%       Color: specifies a MATLAB-defined color, or a 3 element vector
%              specifying an RGB value to fill in the polyhedron faces
%       EdgeColor: specifies a MATLAB-defined color, or a 3 element vector
%              specifying an RBG value to color the polyhedron lines
%       EdgeThickness: Line thickness of the edges
arguments
    F
    V
    Color=[211/255 211/255 211/255];
    EdgeColor=[0 0 0];
    EdgeThickness = 0.5;
end

patch('Faces',F,'Vertices',V,...
    'LineWidth',EdgeThickness,'EdgeColor',EdgeColor,'FaceColor',Color);
view(3)
axis equal

end


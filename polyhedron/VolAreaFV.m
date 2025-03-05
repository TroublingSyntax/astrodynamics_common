function [Volume,Area] = VolAreaFV(faces,vertices)

% <strong>VolAreaFV</strong> - Volume and area for a given triangulated shape.
%
%   [Volume,Area] = VolAreaFV(faces,vertices)
%
% -------------------------------------------------------------------------
%   <strong>INPUTS      unit        format          Description</strong>
% -------------------------------------------------------------------------
%   faces       [-]         (nfaces x 3)    vertices id for each face (row)
%   vertices    [L]         (nvertices x 3) coordinates of each vertex
%                                           (row) in the unit provided by
%                                           the input file
%
% -------------------------------------------------------------------------
%   <strong>OUTPUTS      unit       format          Description</strong>
% -------------------------------------------------------------------------
%   Volume      [L^3]       scalar          volume in the unit provided by
%                                           the input file
%   Area        [L^2]       scalar          area in the unit provided by
%                                           the input file
% -------------------------------------------------------------------------
%
%   The function computes volume enclosed within a triangulated shape and
%   its surface area.
%
% ----<strong>EXAMPLE</strong>--------------------------------------------------------------
%       % define inputs
%       FV=shapefile2FV('asteroid.txt','vf_txt',1);
%       faces=FV.faces;
%       vertices=FV.vertices;
%       % call the function
%       [Volume,Area] = VolAreaFV(faces,vertices);
% -------------------------------------------------------------------------
%
% ###############################################
% #                                             #
% #     <strong>AstroPoliMi MatLab Library [TBC]</strong>        #
% #                                             #
% #     Politecnico di Milano                   #
% #     Aerospace Science & Technology Dept.    #
% #     Via La Masa 34, 20156 Milano-Italia     #
% #     http://www.websiteaddress.[TBD]         #
% #     <a href="mailto:contact.email@domain.TBD">contact.email@domain.[TBD]</a>              #
% #                                             #
% ###############################################
%
% Function :  VolAreaFV (Version 1.0)
% Author(s):  Fabio Ferrari (<a href="mailto:fabio1.ferrari@polimi.it">fabio1.ferrari@polimi.it</a>)
% 
% <strong>Changelog</strong>:
%       2015-10-28 - v1.0 - Fabio Ferrari
%                     - First release.
%
% This code comes with no guarantee or warranty of any kind. [PROVISIONAL]
%
% [Licence policies, TBD]


% edges 3->2
e32 = vertices(faces(:,2),1:3) - vertices(faces(:,3),1:3);
% edges 2->1
e21 = vertices(faces(:,1),1:3) - vertices(faces(:,2),1:3);

% face normals (non-unit vector)
Fn = cross(e32, e21, 2);

% area of each triangle
areaTriang = 0.5*sqrt(Fn(:,1).^2+Fn(:,2).^2+Fn(:,3).^2);
% total area
Area = sum(areaTriang);

% face normals (unit vector)
Fn = bsxfun(@rdivide, Fn, sqrt(sum(Fn.^2, 2)));    
% z component of normal unit vector for each triangle
nz = -Fn(:,3);
% mean of z components of vertices
zMean = (vertices(faces(:,1),3) + vertices(faces(:,2),3) + vertices(faces(:,3),3)) / 3;

% volume of each triangle
volTriang = areaTriang.*zMean.*nz;
% total volume (divergence theorem)
Volume = sum(volTriang);


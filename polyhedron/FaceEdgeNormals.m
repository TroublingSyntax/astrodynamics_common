function [Fn,En] = FaceEdgeNormals(faces,vertices)

% <strong>FaceEdgeNormals</strong> - Normal unit vectors of faces and edges in a 3D mesh.
%
%   [Fn,En] = FaceEdgeNormals(faces,vertices)
%
% -------------------------------------------------------------------------
%   <strong>INPUTS      unit        format          Description</strong>
% -------------------------------------------------------------------------
%   faces       [-]         (nfaces x 3)    vertices id for each face (row)
%   vertices    [L]         (nvertices x 3) coordinates of each vertex
%                                           (row) in the unit provided
%
% -------------------------------------------------------------------------
%   <strong>OUTPUTS      unit       format          Description</strong>
% -------------------------------------------------------------------------
%   Fn          [-]         (nfaces x 3)    normal vector to each face 
%                                           (row)
%   En          [-]         (nfaces x 9)    normal vector to each edge of 
%                                           each face (row)
%                                           [e12 e23 e31]
%                                           [    ...    ]
% -------------------------------------------------------------------------
%
%   The function computes all vector normals of faces and edges of a 3D
%   mesh given in faces-vertices format.
%   The orientation of the normal vector to face A is defined by the sign 
%   of cross product between edges 1->2 and 1->3, with 1, 2 and 3 being 
%   vertices of A.
%   The orientation of the normal vector to edge 1->2 in face A is defined 
%   by the sign of cross product between edge 1->2 and normal to face A
%   (couterclockwise).
%   See (1) for details on notation.
%
% ----<strong>EXAMPLE</strong>--------------------------------------------------------------
%       % define inputs
%       FV = shapefile2FV('asteroid.txt','vf_txt',1);
%       % call the function
%       [Fn,En] = FaceEdgeNormals(FV.faces,FV.vertices);
%       % plot the results
%       TR = triangulation(FV.faces,FV.vertices);
%       CNT=incenter(TR);
%       hold on
%       quiver3(CNT(:,1),CNT(:,2),CNT(:,3),Fn(:,1),Fn(:,2),Fn(:,3))
% -------------------------------------------------------------------------
% 
% <strong>References</strong>
%   (1) Werner, RA and Scheeres, DJ (1997) Exterior Gravitation of a
%   Polyhedron derived and compared with Harmonic and Mascon Gravitation
%   representations of Asteroid 4769 Castalia.
%   Celestial Mechanics and Dynamical Astronomy 65: 313-344.
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
% Function :  FaceEdgeNormals (Version 1.0)
% Author(s):  Fabio Ferrari (<a href="mailto:fabio1.ferrari@polimi.it">fabio1.ferrari@polimi.it</a>)
% 
% <strong>Changelog</strong>:
%       2015-10-23 - v1.0 - Fabio Ferrari
%                     - First release.
%
% This code comes with no guarantee or warranty of any kind. [PROVISIONAL]
%
% [Licence policies, TBD]



% edges 1->2
    e12 = vertices(faces(:,2),1:3) - vertices(faces(:,1),1:3);
    % normalization
    e12 = bsxfun(@rdivide, e12, sqrt(sum(e12.^2, 2)));
% edges 1->3
    e13 = vertices(faces(:,3),1:3) - vertices(faces(:,1),1:3);
    % normalization
    e13 = bsxfun(@rdivide, e13, sqrt(sum(e13.^2, 2)));
% edges 2->3
    e23 = vertices(faces(:,3),1:3) - vertices(faces(:,2),1:3);
    % normalization
    e23 = bsxfun(@rdivide, e23, sqrt(sum(e23.^2, 2)));
    
% face normals
    Fn = cross(e12, e13, 2);
    % normalization
    Fn = bsxfun(@rdivide, Fn, sqrt(sum(Fn.^2, 2)));

% edge normals
    % normals to edges 1->2
        En12 = cross(e12, Fn, 2);
        % normalization
        En12 = bsxfun(@rdivide, En12, sqrt(sum(En12.^2, 2)));
    % normals to edges 2->3
        En23 = cross(e23, Fn, 2);
        % normalization
        En23 = bsxfun(@rdivide, En23, sqrt(sum(En23.^2, 2)));
    % normals to edges 3->1
        En31 = cross(-e13, Fn, 2);
        % normalization
        En31 = bsxfun(@rdivide, En31, sqrt(sum(En31.^2, 2)));
    % assembling En
    En = [En12 En23 En31];

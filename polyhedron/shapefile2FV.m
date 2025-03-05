function [FV] = shapefile2FV(filename,filetype,plotflag)

% <strong>shapefile2FV</strong> - Conversion from shape input file to faces and vertices.
%
%   [FV] = shapefile2FV(filename,filetype,plotflag)
%
% -------------------------------------------------------------------------
%   <strong>INPUTS      unit        format          Description</strong>
% -------------------------------------------------------------------------
%   filename    [-]         string          name of the file
%   filetype    [-]         string          type of the file ('vf_txt' or 
%                                           'tab')
%   plotflag    [-]         scalar          boolean: 1 or 0 to show or not 
%                                           3D plot
%
% -------------------------------------------------------------------------
%   <strong>OUTPUTS      unit       format          Description</strong>
% -------------------------------------------------------------------------
%   FV          [-]         structure       structure with fields:
%   FV.faces    [-]         (nfaces x 3)    vertices id for each face (row)
%   FV.vertices [L]         (nvertices x 3) coordinates of each vertex
%                                           (row) in the unit provided by
%                                           the input file
% -------------------------------------------------------------------------
%
%   The function converts the shape input file into faces-vertices format.
%   List of filetypes
%       'vf_txt'
%       ------------------------
%      |    v     X   Y   Z     |
%      |    ...                 |
%      |    f     vid vid vid   |
%      |    ...                 |
%       ------------------------
%       'tab'
%       ------------------------
%      |    1     X   Y   Z     |
%      |    ...                 |
%      |    1     vid vid vid   |
%      |    ...                 |
%       ------------------------
%
% ----<strong>EXAMPLE</strong>--------------------------------------------------------------
% [Minimal working example of function use, with view of the results.]
%       % define inputs
%       filename='asteroid.txt';
%       filetype='vf_txt';
%       plotflag=1;
%       % call the function
%       [FV] = shapefile2FV(filename,filetype,plotflag);
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
% Function :  shapefile2FV (Version 1.0)
% Author(s):  Fabio Ferrari (<a href="mailto:fabio1.ferrari@polimi.it">fabio1.ferrari@polimi.it</a>)
% 
% <strong>Changelog</strong>:
%       2015-10-23 - v1.0 - Fabio Ferrari
%                     - First release.
%
% This code comes with no guarantee or warranty of any kind. [PROVISIONAL]
%
% [Licence policies, TBD]



switch filetype
    
    case 'vf_txt'
        file_id=fopen(filename);
        Cdata=textscan(file_id,'%c %f %f %f');
        fclose(file_id);
        FVdata=[Cdata{2} Cdata{3} Cdata{4}];
        FV.vertices=FVdata(Cdata{1}=='v',:);
        FV.faces=FVdata(Cdata{1}=='f',:);
        if min(min(FV.faces))==0
            FV.faces = FV.faces + 1;
        end
        
    case 'tab'
        file_id=fopen(filename);
        Cdata=textscan(file_id,'%d %f %f %f');
        fclose(file_id);
        nvertex=Cdata{1}(1);
        FV.vertices=[Cdata{2}(2:1+nvertex) Cdata{3}(2:1+nvertex) Cdata{4}(2:1+nvertex)];
        FV.faces=[Cdata{2}(2+nvertex:end) Cdata{3}(2+nvertex:end) Cdata{4}(2+nvertex:end)];

end



% plot the mesh
if plotflag
    
    ast=patch('Vertices',FV.vertices,'Faces',FV.faces);
    set(ast,'Edgealpha',0.5);
    set(ast,'EdgeColor','white');
    lighting gouraud
    axis equal
    
end

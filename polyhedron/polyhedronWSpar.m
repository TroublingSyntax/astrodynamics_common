function [U,DU,DDU,D2U] = polyhedronWSpar(FP,faces,vertices,density,plotflag)

% <strong>polyhedronWSpar</strong> - Gravity effects of polyhedron model (parallel computing 
%                                    version)
%
%   [U,DU,DDU,D2U] = polyhedronWSpar(FP,faces,vertices,density,plotflag)
%
% -------------------------------------------------------------------------
%   <strong>INPUTS      unit        format          Description</strong>
% -------------------------------------------------------------------------
%   FP          [m]         (3x1) vector    field point
%   faces       [-]         (nfaces x 3)    vertices id for each face (row)
%   vertices    [m]         (nvertices x 3) coordinates of each vertex
%                                           (row)
%   density     [kg/m^3]    scalar          constant density of body
%   plotflag    [-]         scalar          boolean: 1 or 0 to show or not 
%                                           3D plot
%
% -------------------------------------------------------------------------
%   <strong>OUTPUTS      unit       format          Description</strong>
% -------------------------------------------------------------------------
%   U           [m^2/s^2]   scalar          gravitational potential
%   DU          [m/s^2]     (3x1) vector    gradient of potential
%                                           (acceleration)
%   DDU         [s^-2]      (3x3) matrix    gradient of gradient of 
%                                           potential
%   D2U         [s^-2]      scalar          laplacian of potential
% -------------------------------------------------------------------------
%
%   The function computes the exact gravitational effects due to a constant
%   density polyhedron.
%   The gravity effects are computed by including the contributions due to
%   edges and faces of the polyhedral model of the body.
%   The code is valid for triangulated mesh in input.
%
% ----<strong>EXAMPLE</strong>--------------------------------------------------------------
%       % define inputs
%       FP=rand(3,1)*1e3;
%       FV=shapefile2FV('asteroid.txt','vf_txt',0);
%       density=2e3;
%       % initialize parallel pool
%       delete(gcp('nocreate'));
%       mypool=parpool(8);
%       % call the function
%       [U,DU,DDU,D2U] = polyhedronWSpar(FP,FV.faces,FV.vertices,density,1);
%       % delete parallel pool
%       delete(mypool)
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
% #     <strong>ASTRA Team MatLab Library</strong>               #
% #                                             #
% #     Politecnico di Milano                   #
% #     Aerospace Science & Technology Dept.    #
% #     Via La Masa 34, 20156 Milano-Italia     #
% #     http://www.websiteaddress.[TBD]         #
% #     <a href="mailto:contact.email@domain.TBD">contact.email@domain.[TBD]</a>              #
% #                                             #
% ###############################################
%
% Function :  polyhedronWSpar (Version 1.0)
% Author(s):  Fabio Ferrari (<a href="mailto:fabio1.ferrari@polimi.it">fabio1.ferrari@polimi.it</a>)
% 
% <strong>Changelog</strong>:
%       2015-10-28 - v1.0 - Fabio Ferrari
%                     - First release.
%
% This code comes with no guarantee or warranty of any kind. [PROVISIONAL]
%
% [Licence policies, TBD]


% general data
G_univ=6.67384*10^-11;      %[m^3 kg^-1 s^-2]


% field point in the body fixed frame
% transpose FP into column vector if not already
[rP,~]=size(FP);
if rP==1
    FP=FP';
end
clear rP

% faces and edges normals
[Fn,En] = FaceEdgeNormals(faces,vertices);

% edges IDs
Edg = edges (triangulation(faces,vertices));


% for loop (each edge)

% number of edges
[nedges,~]=size(Edg);
% initialization of edges terms
Uj_edge=zeros(1,nedges);
DUj_edge=zeros(3,nedges);
DDUj_edge=zeros(3,3,nedges);

parfor w=1:nedges  % for each edge
    
    % vertices ids
    v_id_i=Edg(w,1);
    v_id_j=Edg(w,2);

    % identification of face A and B
    [Fi,~]=find(faces==v_id_i);
    [Fj,~]=find(faces==v_id_j);
    q=1;
    FAind=[];
    FBind=[];
    while isempty(FBind)
        if isempty(FAind)
            FAind=find(Fi==Fj(q));
        else
            FBind=find(Fi==Fj(q));
        end
        q=q+1;
    end
    FA=Fi(FAind);
    FB=Fi(FBind);
    
    % normal to face A
    nA=Fn(FA,:)';
    nA=nA/norm(nA);

    % normal to face B
    nB=Fn(FB,:)';
    nB=nB/norm(nB);

    % find edge ij in face A
    i1A=find(faces(FA,:)==v_id_i);
    i2A=find(faces(FA,:)==v_id_j);
    if (i1A==1 && i2A==2) || (i1A==2 && i2A==1)
        spanA=1:3;
    elseif (i1A==2 && i2A==3) || (i1A==3 && i2A==2)
        spanA=4:6;
    else
        spanA=7:9;
    end
    % normal to edge (outgoing from face A)
    nijA=En(FA,spanA)';
    nijA=nijA/norm(nijA);

    % find edge ji in face B
    i1B=find(faces(FB,:)==v_id_i);
    i2B=find(faces(FB,:)==v_id_j);
    if (i1B==1 && i2B==2) || (i1B==2 && i2B==1)
        spanB=1:3;
    elseif (i1B==2 && i2B==3) || (i1B==3 && i2B==2)
        spanB=4:6;
    else
        spanB=7:9;
    end
    % normal to edge ji (outgoing from face B)
    njiB=En(FB,spanB)';
    njiB=njiB/norm(njiB);

    % coordinates of vertices i and j
    Pi=vertices(v_id_i,:)';
    Pj=vertices(v_id_j,:)';
    % length of the edge ij
    eij=norm(Pj-Pi);
    % vectors Ri and Rj and their norms ri and rj
    Ri=Pi-FP;
    Rj=Pj-FP;
    ri=norm(Ri);
    rj=norm(Rj);

    % term Le
    Le = log( (ri+rj+eij)/(ri+rj-eij) );
    % term Ee
    Ee = nA*nijA' + nB*njiB';
    % term re
    re = Ri;

    % edge term in the potential function
    Uj_edge(w) = re'*Ee*re * Le;

    % edge term in the gradient of potential function
    DUj_edge(:,w) = Ee*re * Le;

    % edge term in the gradient of the gradient of potential function
    DDUj_edge(:,:,w) = Ee * Le;

end

% sum over the edges
U_edge  = sum(Uj_edge);
DU_edge = sum(DUj_edge,2);
DDU_edge= sum(DDUj_edge,3);


% for loop (each face)

% number of faces
[nfaces,~]=size(faces);
% initialization of potential function
Uj_face=zeros(1,nfaces);
DUj_face=zeros(3,nfaces);
DDUj_face=zeros(3,3,nfaces);
D2Uj_face=zeros(1,nfaces);

parfor h=1:nfaces  % for each face
    
    % initialization (vertices 1 2 3 corresponds to columns 1 2 3)
    P123=zeros(3,3);
    R123=zeros(3,3);
    for k=1:3   % vertices
        % vertices id of the hth triangular face
        v_id=faces(h,k);
        % vertices coordinates in the body fixed frame (column vectors)
        P123(:,k)=vertices(v_id,:)';
        % vectors from field point to vertices
        R123(:,k)=P123(:,k)-FP;
    end
    
    % normal to face
    nF=Fn(h,:)';
    nF=nF/norm(nF);
    
    % vectors R1, R2, R3
    R1=R123(:,1);
    R2=R123(:,2);
    R3=R123(:,3);
    r1=norm(R1);
    r2=norm(R2);
    r3=norm(R3);
    
    % term omegaf
    NUMomegaf = R1'*cross(R2,R3);
    DENomegaf = r1*r2*r3 + r1*(R2'*R3) + r2*(R3'*R1) + r3*(R1'*R2);
    omegaf = 2*atan(NUMomegaf/DENomegaf);
    % term Ff
    Ff = nF*nF';
    % term rf
    rf = R1;
    
    % potential of a single triangular face
    Uj_face(h) = rf'*Ff*rf * omegaf;
    
    % gradient of potential of a single triangular face
    DUj_face(:,h) = Ff*rf * omegaf;
    
    % gradient of gradient of potential of a single triangular face
    DDUj_face(:,:,h) = Ff * omegaf;
    
    % laplacian of potential of a single triangular face
    D2Uj_face(h) = omegaf;

end

% sum over the faces
U_face   = sum(Uj_face);
DU_face  = sum(DUj_face,2);
DDU_face = sum(DDUj_face,3);
D2U_face = sum(D2Uj_face);


% put all components together

% potential
U = 0.5*density*G_univ * ( U_edge - U_face );

% gradient of potential
DU = density*G_univ *  ( - DU_edge + DU_face );

% gradient of gradient of potential
DDU = density*G_univ * ( DDU_edge - DDU_face );

% laplacian of potential
D2U = - density*G_univ * D2U_face ;


% plot the mesh
if plotflag
    ast=patch('Vertices',vertices,'Faces',faces);
    set(ast,'Edgealpha',0.5);
    set(ast,'EdgeColor','white');
    % set(ast,'facealpha',0.2);
    % set(ast,'FaceColor','white');
    lighting gouraud
    axis equal
end


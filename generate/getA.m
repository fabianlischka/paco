% stokes.edp provides: u.txt, v.txt, mass.txt, rih.txt, stiff.txt,
%                      and the mesh: stokes.msh
% where u, v 1370 x 1; mass, stiff 736 x 736; stokes.msh in mesh format
% inputs used here appear to be: 
% ff, vv; u, v; mi, mj, ms (=mass); ki, kj, ks (=?); bj, bi, ones (=?)
% nu
% where:
% u, v = 1370 x 1   (from u, v)
% ff   = 1370 x 3   (from mesh), int v1, v2, v3: indices into
%                    vv (vertices))
% vv   =  736 x 3   (from mesh), float x, y; int l_r: coords and type of
%                    vertices)
% mi,j,s = 4946 x 1 (from mass.txt)
% ki,j,s = 4946 x 1 (from stiff.txt)
% bi,j,s =   41 x 1 (from Rih.txt)

% [ik,jk,sk] = fvm( ff(k,:), vv(ff(k,:),1:2), [u(k), v(k)], vv(ff(k,:),3) );
%                   -> vv    -> q             -> u          -> where

ii = [];
jj = [];
ss = [];



% vv = vertices: coordx, coordy, type
% ff = faces: vertex1, vertex2, vertex3
% geometry, comes from stokes.msh



load bay;
nu = 0.01;

for k=1:length(ff),
    [ik,jk,sk] = fvm(ff(k,:),vv(ff(k,:),1:2),[u(k), v(k)], vv(ff(k,:),3));
    ii = [ii; ik]; % append
    jj = [jj; jk];
    ss = [ss; sk];
end;

A = sparse(ii,jj,ss,size(vv,1),size(vv,1),size(ii,1));   % A
M = sparse(mi,mj,ms,size(vv,1),size(vv,1));              % 
M2 = spdiags(sqrt(diag(M)),0,size(vv,1),size(vv,1));     % 
K = sparse(ki,kj,ks,size(vv,1),size(vv,1));              % 
B = sparse(bj,bi,ones(length(bi),1),size(vv,1),max(bi)); % B

%%
AA = A + nu * K;
dt = 0.5;
Am = M2\AA/M2;  % this is A
Bm = M2\B;      % this is B
Nt = 40;

y = zeros(size(vv,1),Nt+1);

for i=1:Nt,
    if (i > 10 && i < 30),
        rn = y(:,i) + dt*Bm*ones(size(B,2),1);
    else
        rn = y(:,i);
    end;
    y(:,i+1) = (eye(size(A))+dt*Am)\rn;
end;
y = M2\y;
%%
for i=1:Nt,
    disp(i)
    patch('Faces',ff,'Vertices',vv(:,1:2),'FaceColor','interp','FaceVertexCData',y(:,i),'LineStyle','none');
    axis equal
    ax=gca;
    ax.CLim=[0 4];
    colorbar;
    axis off
    pause;
end;

function g = gifti_inflator(surf,infl)
% surf is a gifti surface (or struct w/ vertices & faces)
% infl is amount of inflation between 0-1
% g    is surf but with inflated vertices
%
% AS2017

if infl > 1; infl = 1/infl; end

g  = gifti(surf);
f  = g.faces;
v  = g.vertices;%/( max(g.vertices(:))*(pi/2) );
n  = vsmooth(v, double(f), infl);

g.vertices = n;
g.faces    = double(g.faces);
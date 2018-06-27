function h = meshmesh(g,C,alph)

f = g.faces;
v = g.vertices;
h = patch('faces',f,'vertices',v);

if nargin < 2 || isempty(C)
    C = [.5 .5 .5];
end
if nargin < 3 || isempty(alph)
    alph = .5;
end

set(h,'FaceColor',C);
box off;
grid off; 
set(h,'EdgeColor','none')
set(gca,'visible','off');
%alpha(alph);
set(h,'FaceAlpha',alph);
function meshmesh(g)

f = g.faces;
v = g.vertices;
h = patch('faces',f,'vertices',v);

C = .5;

set(h,'FaceColor',[C C C]);
box off;
grid off; 
set(h,'EdgeColor','none')
alpha(.3);
set(gca,'visible','off');
function meshmesh(g)

f = g.faces;
v = g.vertices;

h = patch('faces',f,'vertices',v);

 set(h,'FaceColor',[.4 .4 .4]);
 box off;
 grid off; 
 set(h,'EdgeColor','none')
 alpha(.2);
 set(gca,'visible','off');
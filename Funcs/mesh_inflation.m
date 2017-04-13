function mesh_inflation(D,i)
% D is inverted spm12 meeg dataset [with meshdata in forward substructure]
% i is between 0 - 1 [def .02]
% 
% AS16

if iscell(D); D = loadarrayspm(D);D=D{1}; end
if nargin == 1; i = .02; end

vert  = D.inv{end}.forward(end).mesh.vert;
face  = D.inv{end}.forward(end).mesh.face;
face  = double(face);

nvert = vsmooth(vert, face, i);
h = patch('faces',face,'vertices',nvert);

set(h,'FaceColor',[.4 .4 .4]);
box off;
grid off;
whitebg(1,'w'); 
camlight('right')
%axis tight
set(h,'EdgeColor','none')
material dull
alpha(.2);
lighting phong
set(gca,'visible','off');
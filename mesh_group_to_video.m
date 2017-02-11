
f = dir('baw*');f={f.name};

for i = 1:length(f); 
    clear D; 
    
    D    = spm_eeg_load(f{i}); % load D
    y{i} = mesh_pca1(D,'Dev'); % smooth until ncomps
end

yy    = cat(3,y{:});
mesh  = D.inv{1}.forward(end).mesh;
mesh  = export(gifti(mesh));

for i = 1:size(yy,2)
    my(:,i,:) = spm_mesh_smooth(mesh,squeeze(yy(:,i,:)),8);
end

mmy = squeeze(mean(my,3));

%sy = HighResMeanFilt(yy,1,4);  % smooth over subs
%my = squeeze(mean(sy,3));

%[C,A,W] = fastica(my);
%W       = pinv(W);

% for i = 1:size(C,1)
%     component{i} = ( (C(i,:)'*W(:,i)') );
% end

time = D.inv{1}.inverse.pst;

plotmeshfov(D,abs(component{1}),5e-4,.3,'Test',[80,45]);
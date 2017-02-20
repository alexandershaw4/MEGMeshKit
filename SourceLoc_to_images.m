

% add tools
addpath(genpath('/home/as08/old_spm12'));
%addpath('/home/as08/MNI2TS/')
%cd /imaging/as08/psp_faces/ ;
cd /Volumes/Extra/Rov_ICA/

%load('datasets')
F = dir('bawin*.mat'); F = {F.name};
%F = strrep(F,'1_30','30_90');
f = loadarrayspm(F);
f = f';

toi = [0 .3];
foi = [];

trans = .8;          % overlay alpha (0 - 1)
typ   = 'trials';

% calculate group responses
% [Neutm,Neut]  = plotmesh_fo_grp_pca(f,[],{'NeutralAll'},toi,foi,typ);
% [Happm,Happ]  = plotmesh_fo_grp_pca(f,[],{'HappyAll'},toi,foi,typ);
% [Angrm,Angr]  = plotmesh_fo_grp_pca(f,[],{'AngryAll'},toi,foi,typ);
% [FTm,FT]      = plotmesh_fo_grp_pca(f,[],{'FTAll'},toi,foi,typ);

[Devm,Dev]   = plotmesh_fo_grp_pca(f,{'Dev'});
[rep1m,rep1] = plotmesh_fo_grp_pca(f,{'rep1'});
[rep2m,rep2] = plotmesh_fo_grp_pca(f,{'rep2'});
[rep3m,rep3] = plotmesh_fo_grp_pca(f,{'rep3'});
[rep4m,rep4] = plotmesh_fo_grp_pca(f,{'rep4'});
[rep5m,rep5] = plotmesh_fo_grp_pca(f,{'rep5'});
[rep6m,rep6] = plotmesh_fo_grp_pca(f,{'rep6'});
[rep7m,rep7] = plotmesh_fo_grp_pca(f,{'rep7'});
[rep8m,rep8] = plotmesh_fo_grp_pca(f,{'rep8'});
[rep9m,rep9] = plotmesh_fo_grp_pca(f,{'rep9'});

cond = {'Dev', 'rep1', 'rep2', 'rep3', 'rep4', 'rep5', 'rep6', 'rep7', 'rep8','rep9'};
% cond = {'NeutralAll' 'HappyAll' 'AngryAll' 'FTAll'};


% write out images for spm
for s = 1:length(f)

    fn = f{s}.fname;
        
    fnam = ['Dev_' fn];
    writeims(f{s},Dev(s,:),fnam);
    
    fnam = ['Rep1_' fn];
    writeims(f{s},rep1(s,:),fnam);
    
    fnam = ['Rep2_' fn];
    writeims(f{s},rep2(s,:),fnam);
    
    fnam = ['Rep3_' fn];
    writeims(f{s},rep3(s,:),fnam);
    
    fnam = ['Rep4_' fn];
    writeims(f{s},rep4(s,:),fnam);
    
    fnam = ['Rep5_' fn];
    writeims(f{s},rep5(s,:),fnam);

    fnam = ['Rep6_' fn];
    writeims(f{s},rep6(s,:),fnam);
    
    fnam = ['Rep7_' fn];
    writeims(f{s},rep7(s,:),fnam);
    
    fnam = ['Rep8_' fn];
    writeims(f{s},rep8(s,:),fnam);
    
    fnam = ['Rep9_' fn];
    writeims(f{s},rep9(s,:),fnam);
    
end


% add tools
addpath(genpath('/home/as08/old_spm12'));
addpath('/home/as08/MNI2TS/')
cd /imaging/as08/psp_faces/ ;

load('datasets')
F = strrep(F,'1_30','30_90');
f = loadarrayspm(F);
f = f';

toi = [0 .35];
foi = [];

trans = .8;          % overlay alpha (0 - 1)
typ   = 'trials';

% calculate group responses
[Neutm,Neut]  = plotmesh_fo_grp_pca(f,[],{'NeutralAll'},toi,foi,typ);
[Happm,Happ]  = plotmesh_fo_grp_pca(f,[],{'HappyAll'},toi,foi,typ);
[Angrm,Angr]  = plotmesh_fo_grp_pca(f,[],{'AngryAll'},toi,foi,typ);
[FTm,FT]      = plotmesh_fo_grp_pca(f,[],{'FTAll'},toi,foi,typ);

cond = {'NeutralAll' 'HappyAll' 'AngryAll' 'FTAll'};

% write out images for spm
for s = 1:length(f)

    fn = f{s}.fname;
        
    fnam = ['Neutral_' fn];
    writeims(f{s},Neut(s,:),fnam);
    
    fnam = ['Happy_' fn];
    writeims(f{s},Happ(s,:),fnam);
    
    fnam = ['Angry_' fn];
    writeims(f{s},Angr(s,:),fnam);
    
    fnam = ['FT_' fn];
    writeims(f{s},FT(s,:),fnam);
        
end
function writeims(D,ol,sname)
% Utility for saving functional mesh overlay to gifti
%
% D  is a subjects' meeg object (for getting templt mesh verts)
% ol is the overlay you want to write [1xnverts]
% sname is filename
%
% [based on spms mesh 2 voxels]
% AS

val          = D.val;
woi          = [0 .3];
foi          = [];
Nw           = size(woi,1);
smooth       = 8;
fmt          = 'mesh'; % or 'img';

sMRIfile = fullfile(spm('dir'),'canonical','avg152T2.nii');
Vin      = spm_vol(sMRIfile);
m        = export(gifti(D.inv{val}.mesh.tess_mni),'patch');
nd       = D.inv{val}.inverse.Nd;
GL       = spm_mesh_smooth(m);
GW       = {ol};
Ne       = ones(1, numel(GW));
Nj       = numel(GW)/Nw;

k  = 1;
iw = [];
ie = [];

for c = 1:length(GW)

        cGW = GW{c};

    for t = 1:Ne(c)
        %-Smooth on the cortical surface
        %------------------------------------------------------------------
        ssq{k} = full(sparse(D.inv{val}.inverse.Is,1,cGW,nd,1));
        ssq{k} = spm_mesh_smooth(GL,ssq{k},smooth);
        
        %-Compute (truncated) moment
        %------------------------------------------------------------------
        lss        = log(ssq{k} + eps);
        i          = lss > (max(lss) - log(32));
        meanlss(k) = mean(lss(i));
        
        iw(k) = c;
        ie(k) = t;
        
        k = k + 1;
    end
end

scale = exp(mean(meanlss));

%-Save as meshes
%==========================================================================
if strcmpi(fmt,'mesh')
    [pth,name] = fileparts(D.fname);
    tag = cell(1,Nw);
    for i = 1:Nw
        tag{i} = ['t' sprintf('%g_', woi(i,:)) 'f' sprintf('%d_', foi)];
    end
    
    %-Save mesh topology
    %----------------------------------------------------------------------
    g = gifti(m);
    save(g,fullfile(pwd,[sname '.surf.gii']));
    
    %-Save mesh data
    %----------------------------------------------------------------------
    for c = 1:numel(ssq)
        
%         fprintf('%s%30s',repmat(sprintf('\b'),1,30),...
%             sprintf('...mesh %d/%d',c,numel(ssq)));                     %-#
        
        %-Filename
        %------------------------------------------------------------------
        con       = mod(iw(c) - 1, Nj) + 1;
        str       = tag{ceil(iw(c)/Nj)};
        bt = '';
        fname     = fullfile(pwd,...
            sprintf('%s_%.0f_%s%.0f%s.gii', sname, val, str, con, bt));
        
        %-Normalise
        %------------------------------------------------------------------
        Contrast = ssq{c} / scale;
        
        %-Write mesh
        %------------------------------------------------------------------
        g = gifti;
        g.cdata = Contrast;
        g.private.metadata(1).name = 'SurfaceID';
        g.private.metadata(1).value = [name '.surf.gii'];
        save(g, fname, 'ExternalFileBinary');
        
        %-Store filename
        %------------------------------------------------------------------
        D.inv{val}.contrast.fname{c} = fname;
    end
    
    %return;
end

%-Normalise and embed in 3D-space
%==========================================================================
tag = cell(1,Nw);
for i = 1:Nw
    tag{i} = ['t' sprintf('%d_', woi(i,:)) 'f' sprintf('%d_', foi)];
end

con       = mod(iw - 1, Nj) + 1;
win       = ceil(iw/Nj);
ucon      = unique(con);
uwin      = unique(win);

n = 0;
for c = 1:length(ucon)
    for w = 1:numel(uwin)
        
        ind = find((con == ucon(c)) & (win == uwin(w)));
        str = tag{uwin(w)};
        
        fname     = fullfile(pwd,...
            sprintf('%s_%.0f_%s%.0f.nii', sname, val, str, ucon(c)));
        
        %-Initialise image
        %----------------------------------------------------------------------
        N      = nifti;
        N.dat  = file_array(fname, [Vin.dim, length(ind)], 'FLOAT32-LE');
        N.mat  = Vin.mat;
        N.mat0 = Vin.mat;
        create(N);
        
        for i = 1:length(ind)
            
            n = n+1;
            
            %-Normalise
            %----------------------------------------------------------------------
            Contrast = ssq{ind(i)} / scale;
            
            %-Interpolate those values into voxels
            %----------------------------------------------------------------------
            RECimage = spm_mesh_to_grid(m, Vin, Contrast);
            
            %-3D smoothing and thresholding
            %----------------------------------------------------------------------
            try spm_smooth(RECimage, RECimage, 1); end
            RECimage = RECimage.*(RECimage > exp(-8));
            
            %-Write (smoothed and scaled) image
            %----------------------------------------------------------------------            
            N.dat(:, :, :, i) = RECimage;  
            
        end
    end
end
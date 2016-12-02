function out = rebuild(D,woi,type,foi)
% Plot source reconstructed MEEG data from spm.
%
%
%
%
% AS

model = D.inv{1};

woi(:,1) = max(woi(:,1),model.inverse.woi(1));
woi(:,2) = min(woi(:,2),model.inverse.woi(2));

if ~isempty(foi); fprintf('using frequency specific responses\n'); end

% inversion parameters
%--------------------------------------------------------------------------
J    = model.inverse.J;                        % Trial average MAP estimate
T    = model.inverse.T;                        % temporal projector
U    = model.inverse.U;                        % spatial  projector[s]
Is   = model.inverse.Is;                       % Indices of ARD vertices
Ic   = model.inverse.Ic;                       % Indices of channels
It   = model.inverse.It;                       % Indices of time bins
pst  = model.inverse.pst;                      % peristimulus time (ms)
Nd   = model.inverse.Nd;                       % number of mesh dipoles
Nb   = size(T,1);                              % number of time bins
Nc   = size(U,1);                              % number of channels
Nw   = size(woi,1);                            % number of contrast windows
Nj   = numel(J);                               % number of conditions

try   scale = model.inverse.scale;               % Trial average MAP estimate
catch scale = 1;
end

% time window
%------------------
w    = 1; % only 1 window per run
fwhm = max(diff(woi(w,:)),8);
t    = exp(-4*log(2)*(pst(:) - mean(woi(w,:))).^2/(fwhm^2));
t    = t/sum(t);

    if ~isempty(foi)
        wt = 2*pi*pst(:)/1000;
        W  = [];
        for f = foi(1):foi(end)
            W = [W sin(f*wt) cos(f*wt)];
        end
        W  = diag(t)*W;
        W  = spm_svd(W,1);
    else
        W  = t(:);
    end

TW     = T'*W;
TTW    = T*TW;

% MAP projector
qC     = model.inverse.qC*trace(TTW'*model.inverse.qV*TTW);
qC     = max(qC,0);

try
    trial = model.inverse.trials;
catch
    trial = D.condlist;
end

for i = 1:Nj
        
        % induced or evoked
        %------------------------------------------------------------------
        iw     = (w - 1)*Nj + i;
        CW{iw} = W;
        
        switch(type)
            
            % energy of conditional mean
            %--------------------------------------------------------------
            case{'evoked'}
              
                JW{iw} = J{i}*TW(:,1);
                GW{iw} = sum((J{i}*TW).^2,2) + qC;
        
        
            out.JW = JW;
            out.GW = GW;
                
            %mean energy over trials
            %--------------------------------------------------------------
            case{'induced'}
                
                JW{iw} = sparse(0);
                JWWJ   = sparse(0);
                
                c = D.indtrial(trial{i}, 'GOOD');
                
                %conditional expectation of contrast (J*W) and its energy
                %----------------------------------------------------------
                Nt    = length(c);
                for j = 1:Nt
                    if ~strcmp(D.modality(1,1), 'Multimodal')
                        
                        %unimodal data
                        %--------------------------------------------------
                        Y  = D(Ic{1},It,c(j))*TTW;
                        Y  = U{1}*Y*scale;
                        
                    else
                        
                        %multimodal data
                        %--------------------------------------------------
                        for k = 1:length(U)
                            Y       = D(Ic{k},It,c(j))*TTW;
                            UY{k,1} = U{k}*Y*scale(k);
                        end
                        Y = spm_cat(UY);
                    end
                    
                    MYW    = model.inverse.M*Y;
                    JW{iw} = JW{iw} + MYW(:,1);
                    JWWJ   = JWWJ   + sum(MYW.^2,2);
                    
                end
                
                
                %conditional expectation of total energy (source space GW)
                %----------------------------------------------------------
                JW{iw} = JW{iw}/Nt;
                GW{iw} = JWWJ/Nt + qC;
                
                out.GW = GW;
                out.JW = JW;
                
            case 'trials'
                
                JW{iw} = {};
                JWWJ   = {};
                
                c = D.indtrial(trial{i}, 'GOOD');
                
                %onditional expectation of contrast (J*W) and its energy
                %----------------------------------------------------------
                Nt    = length(c);
                for j = 1:Nt
                    if ~strcmp(D.modality(1,1), 'Multimodal')
                        
                        %unimodal data
                        %--------------------------------------------------
                        Y     = D(Ic{1},It,c(j))*TTW;
                        Y     = U{1}*Y*scale;
                        
                    else
                        
                        %multimodal data
                        %--------------------------------------------------
                        for k = 1:length(U)
                            Y       = D(Ic{k},It,c(j))*TTW;
                            UY{k,1} = U{k}*Y*scale(k);
                        end
                        Y = spm_cat(UY);
                    end
                    
                    MYW       = model.inverse.M*Y;
                    JW{iw}{j} = MYW(:,1);
                    GW{iw}{j} = sum(MYW.^2,2) + qC;
                end
                out.JW = JW;
                out.GW = GW;
        end
        
end
    

function y = mesh_ica1(D,cond,varargin)
% mesh-based ica [with plots] for spm source reconstructed m/eegs
%
% AS

try Strct2WkspLocal(varargin{:}); end   % unpack configuration

try doplot; catch doplot = 0;       end % plot 
try doimg;  catch doimg  = 0;       end % save images [g/nifti]
try invi;   catch invi   = 1;       end % inversion index (see D.val)
try woi;    catch woi    = [0 .35]; end % window of interest (secs)
try foi;    catch foi    = [];      end % opt. foi [but reduces rank]
try k;      catch k      = 8;       end % smoothing kern size
try neig;   catch neig   = 20;      end % num comps 

coni   = strmatch(cond,D.inv{1}.inverse.trials);

%-----------------------------------------------
inv  = D.inv{invi};
vert = inv.forward(end).mesh.vert;
face = inv.forward(end).mesh.face;
pst  = inv.inverse.pst;
T    = inv.inverse.T;
Ja   = inv.inverse.J;
J    = Ja{coni};       
time = D.inv{1}.inverse.pst;
mesh = inv.forward(end).mesh;

% get temp projr for woi
w    = 1; 
fwhm = max(diff(woi(w,:)),8);
t    = exp(-4*log(2)*(pst(:) - mean(woi(w,:))).^2/(fwhm^2));
t    = t/sum(t);
W    = t(:);

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

TW   = T'*W;
TTW  = T*TW;

% full pst:
if isempty(foi); wtime = J*T';
else             wtime = J*(TW*TTW');
end

wi    = [findthenearest(time,woi(1)*1000):findthenearest(time,woi(2)*1000)];
TOI   = abs(wtime(:,wi));

[C, A, W] = fastica(TOI,'numOfIC',neig);
%y     = reduce_eig_mesh(mesh,TOI,neig,k);

nc = size(C,1);
W  = pinv(W);


for i = 1:nc
    this{i} = ( C(i,:)'*W(:,i)' )';   % spatio-temporal component projections
end

for i = 1:nc
    y{i} = this{i}*TTW;           % temporal projector TTW
end

%y = this;

if doplot
    nv = 2;
    nh = ceil(nc/nv);

    for i = 1:nc
        subplot(nv,nh,i),plotmeshfo(D,abs(y{i}),50,.3); 
    end
end


if doimg
    for i = 1:nc
        writeims(D,mean(y{i},2),['comp_' num2str(i) '_' cond]);
    end
end

y.y = y;
y.m = this;

end

function Strct2WkspLocal(S)

f = fieldnames(S);

for i = 1:length(f)
    assignin('caller',f{i},(S.(f{i})));
end


end

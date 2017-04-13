function y = reduce_eig_mesh(mesh,x,d,varargin)
% mesh is a patch structure or gifti
% x is a double matrix
% d is number of comps 
% optional k is smooth kern size
% y is smoothed iteratively until 90% of its variance is explained by d 
% principal components
%
% AS

try k = varargin{1}; catch k = 8; end

[~, n] = PEig90(x);
ncyc   = 0;
maxv   = max(x(:));
en     = 100;

while n > d
    ncyc = ncyc + 1;
    x = spm_mesh_smooth(export(gifti(mesh)), x, k);
    [~, n] = PEig90(x);
    
    en = 90*(en/100); fprintf('\b');fprintf('(%d %%)\n',round(en));
    
    y = (x./max(x(:)))*maxv;
    %x = y;
    
    % Put some limit on cycles
    dn(ncyc) = n;
    if ncyc > 4  && all(dn==n); k = k * 4; end
    if ncyc > 20 && all(dn(end-6:end)==dn(end)); 
       fprintf('Didnt converge: giving up\n'); y = x; return; 
    end
end

fprintf('Finished after %d cycles\n',ncyc);
y = x;

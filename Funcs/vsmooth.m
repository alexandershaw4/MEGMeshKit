function V = vsmooth(Vertices, Faces, VoxSize)
% vertex smoothing

DampingFactor = 0.91;
MaxNormDispl  = 0.5 * VoxSize;
Verbose       = false;
nV            = size(Vertices, 1);
IterTol       = 100;
DisplTol      = 0.01 * VoxSize;
DisplTol      = [DisplTol, 0];


%C = sparse(Faces(:), [Faces(:, 2); Faces(:, 3); Faces(:, 1)], true);
C = Faces';
CCell = cell(nV, 1);
for v = 1:nV
    CCell{v} = find(C(:, v));
end
clear C

V = Vertices;
W = V;
LastMaxDispl = [inf, inf];
Iter = 0;
NormDispl = zeros(nV, 1);
while LastMaxDispl(1) > DisplTol(1) && LastMaxDispl(2) > DisplTol(2) && ...
        Iter < IterTol
    Iter      = Iter + 1;
    [N, VdA]  = VertexNormals(V);
    VWeighted = VdA * [1, 1, 1] .* V;
    
    for v = 1:nV
        NeighdA = sum(VdA(CCell{v}));
        a = sum(VWeighted(CCell{v}, :) / NeighdA, 1); % / nC(v);
        d = (a - V(v, :)) * N(v, :)' / (NeighdA/VdA(v) + 1); % / (nC(v) + 1);
        W(v, :) = W(v, :) + DampingFactor * ( a - V(v, :) - d * N(v, :) );
        W(CCell{v}, :) = W(CCell{v}, :) - DampingFactor * d * N(CCell{v}, :);
    end
    D = NormDispl + dot((W - V), N, 2);
    NormDispl = sign(D) .* min(abs(D), MaxNormDispl);
    D = D - NormDispl;
    Where = abs(D) > DisplTol(1) * 1e-6; % > 0, but ignore precision errors.
    if any(Where)
        W(Where, :) = W(Where, :) - [D(Where), D(Where), D(Where)] .* N(Where, :);
    end
    
    LastMaxDispl(1) = sqrt( max(sum((W - V).^2, 2)) );
    LastMaxDispl(2) = max(dot(W - V, N, 2));
    V = W;
end

if Iter >= IterTol
    warning('SurfaceSmooth did not converge within %d iterations. \nLast max point displacement = %f', ...
        IterTol, LastMaxDispl(1));
elseif Verbose
    fprintf('SurfaceSmooth converged in %d iterations. \nLast max point displacement = %f\n', ...
        Iter, LastMaxDispl(1));
end
if Verbose && IterTol > 0
    [~, ~, FNdA, FdA] = VertexNormals(V);
    FaceCentroidZ = ( V(Faces(:, 1), 3) + ...
        V(Faces(:, 2), 3) + V(Faces(:, 3), 3) ) /3;
    Post.Volume = FaceCentroidZ' * FNdA(:, 3);
    Post.Area = sum(FdA);
    fprintf('Total enclosed volume after smoothing: %g\n', Post.Volume);
    fprintf('Relative volume change: %g %%\n', ...
        100 * (Post.Volume - Pre.Volume)/Pre.Volume);
    fprintf('Total area after smoothing: %g\n', Post.Area);
    fprintf('Relative area change: %g %%\n', ...
        100 * (Post.Area - Pre.Area)/Pre.Area);
end




% ----------------------------------------------------------------------
% Calculate dA normal vectors to each vertex.
    function [N, VdA, FNdA, FdA] = VertexNormals(V)
    N = zeros(nV, 3);
    % Get face normal vectors with length the size of the face area.
    FNdA = CrossProduct( (V(Faces(:, 2), :) - V(Faces(:, 1), :)), ...
        (V(Faces(:, 3), :) - V(Faces(:, 2), :)) ) / 2;
    % For vertex normals, add adjacent face normals, then normalize.  Also
    % add 1/3 of each adjacent area element for vertex area.
    FdA = sqrt(FNdA(:,1).^2 + FNdA(:,2).^2 + FNdA(:,3).^2);
    VdA = zeros(nV, 1);
    for ff = 1:size(Faces, 1) % (This is slow.)
        N(Faces(ff, :), :) = N(Faces(ff, :), :) + FNdA([ff, ff, ff], :);
        VdA(Faces(ff, :), :) = VdA(Faces(ff, :), :) + FdA(ff)/3;
    end
    N = bsxfun(@rdivide, N, sqrt(N(:,1).^2 + N(:,2).^2 + N(:,3).^2));
    end

end

% Much faster than using the Matlab version.
function c = CrossProduct(a, b)
c = [a(:,2).*b(:,3)-a(:,3).*b(:,2), ...
    a(:,3).*b(:,1)-a(:,1).*b(:,3), ...
    a(:,1).*b(:,2)-a(:,2).*b(:,1)];
end









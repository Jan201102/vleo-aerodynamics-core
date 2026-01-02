% Tetrahedron geometry data for MATLAB

% Define tetrahedron vertex positions
v0 = [0.0, -0.5, 0.0];   % Bottom left
v1 = [0.0, 0.5, 0.0];    % Bottom right
v2 = [0.0, 0.0, 0.5];    % Top
v3 = [-0.5, 0.0, 0.0];   % Back

% Tetrahedron vertices - 4 triangular faces (12 vertices total)
% Hardcoded 3x3x4 array: rows = [x;y;z], columns = vertices, pages = faces
V = zeros(3,3,4);

V(:,:,1) = [v0(:), v1(:), v3(:)];   % Face 1: v0, v1, v3
V(:,:,2) = [v1(:), v2(:), v3(:)];   % Face 2: v1, v2, v3
V(:,:,3) = [v2(:), v0(:), v3(:)];   % Face 3: v2, v0, v3
V(:,:,4) = [v0(:), v2(:), v1(:)];   % Face 4: v0, v2, v1

%compute normals
% V is 3x3xF
F = size(V,3);
faceNormals = zeros(3,F);

% V is 3x3xF
F = size(V,3);
centroids = squeeze(mean(V,2));   % results in 3xF


for f = 1:F
    verts = V(:,:,f);        % 3x3: cols are v1,v2,v3
    v1 = verts(:,1);
    v2 = verts(:,2);
    v3 = verts(:,3);
    n = cross(v3 - v1,v2 - v1);   % unnormalized normal (right-hand rule)
    nrm = norm(n);
    if nrm > 0
        faceNormals(:,f) = n / nrm; % unit normal
    else
        faceNormals(:,f) = [0;0;0]; % degenerate face
    end
end


determineShadowedTriangles(V,centroids,faceNormals,wind_dir)
determineShadowedTrianglesGPU(V,faceNormals,wind_dir)


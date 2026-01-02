function ind_shadowed = determineShadowedTrianglesGPU(vertices, normals, dir)
% determineShadowedTriangles - Determine which triangles are shadowed by others
%
%   ind_shadowed = determineShadowedTriangles(vertices, centroids, normals, dir)
%   calculates which triangles are shadowed by at least one of the other triangles
%   along a given a direction dir.
%   A triangle is marked as shadowed if it is not visible in the wind
%   direction
%
%   Inputs:
%   vertices: 3x3xN array of vertices of N triangles, each 3x3 matrix represents on triangle
%             of which each column represents the x, y, z coordinates of one of the 
%             triangle's vertices
%   centroids: 3xN array of surface centroids of N triangles
%   normals: 3xN array of surface normals of N triangles
%   dir: 3x1 array representing the a direction along which the shadowing is determined
%
%   Outputs:
%   ind_shadowed: 1xN logical array indicating which triangles are shadowed
%

%% Principle
% first rear facing triangles are determined, then the visible triangles
% are computed on gpu. lastly both sets (visible and rear-facing) are
% joined and returned

% Angles between flow and normals
delta = real(acos(-dir' * normals));
num_triangles = size(normals,2);
num_vertices = num_triangles*3;

% Determine which triangles are flow-facing and which are rear-facing
ind_flow_facing = (delta < pi/2);
ind_rear_facing = ~ind_flow_facing;
gpu_implementation_path = "Dependencies/shader/BinaryShader/matlab/BinaryShader";
if ~isFolderOnPath(gpu_implementation_path)
    addpath(gpu_implementation_path);
end
if any(ind_flow_facing)
    %reshape vertices to row vector
    shaded = false(1,num_triangles);
    vertices_flat = reshape(vertices, 1, []);
    triangle_ids = zeros(1,num_vertices);
    for id = 1:num_triangles
        triangle_ids(id*3-2:id*3) = id;
    end

    shadedArg = clibConvertArray("clib.BinaryShader.Bool",shaded);
    verticesArg = clibConvertArray("clib.BinaryShader.Float",vertices_flat);
    triangleIDsArg = clibConvertArray("clib.BinaryShader.UnsignedInt",triangle_ids);

    % 4. Call the function
    try
        % Pass the clib objects and the uint64 lengths
        clib.BinaryShader.BinaryRenderer(verticesArg,triangleIDsArg,shadedArg, -dir(1),-dir(2),-dir(3));
    
        % 5. Convert back to MATLAB to see the result
        ind_gpu_visible = logical(shadedArg);
    catch ME
        fprintf('Error: %s\n', ME.message);
    end
    
    ind_visible = ind_rear_facing | ind_gpu_visible;
    ind_shadowed = ~ind_visible;

end
end

function onPath = isFolderOnPath(folder)
% isFolderOnPath  Return true if folder is on the MATLAB search path.
%   folder can be absolute or relative (relative resolved with pwd).

% Resolve to absolute canonical form (remove trailing filesep)
folderAbs = char(java.io.File(folder).getCanonicalPath());  % returns absolute
% Get current MATLAB path entries and canonicalize each
entries = strsplit(path, pathsep);
for k = 1:numel(entries)
    try
        entries{k} = char(java.io.File(entries{k}).getCanonicalPath());
    catch
        % ignore invalid entries
        entries{k} = '';
    end
end
onPath = any(strcmpi(folderAbs, entries));
end

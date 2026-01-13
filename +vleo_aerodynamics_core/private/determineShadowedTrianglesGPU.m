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
%   normals: 3xN array of surface normals of N triangles
%   dir: 3x1 array representing the a direction along which the shadowing is determined
%
%   Outputs:
%   ind_shadowed: 1xN logical array indicating which triangles are shadowed
%

% Declare extrinsic functions for code generation
coder.extrinsic('addpath');
coder.extrinsic('clibConvertArray');
coder.extrinsic('clib.BinaryShader.BinaryRenderer');
coder.extrinsic('logical');

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

if any(ind_flow_facing)
    %normalize geometry to -1,1 in all dimensions
    %1. center geometry
    vertices_centered = vertices - mean(vertices, [2,3]);
    %2. scale geometry
    max_extent = max(abs(vertices_centered), [], 'all');
    vertices_normalized = vertices_centered / max_extent;

    %reshape vertices to row vector
    shaded = false(1,num_triangles);
    vertices_flat = single(reshape(vertices_normalized, 1, []));  % Convert to single precision
    triangle_ids = uint32(zeros(1,num_vertices));  % Use uint32 type
    for id = 1:num_triangles
        triangle_ids(id*3-2:id*3) = id;
    end

    shadedArg = clibConvertArray("clib.BinaryShader.Bool",shaded);
    verticesArg = clibConvertArray("clib.BinaryShader.Float",vertices_flat);
    triangleIDsArg = clibConvertArray("clib.BinaryShader.UnsignedInt",triangle_ids);

    % 4. Call the function
    % Initialize output variable with correct type for code generation
    ind_gpu_visible = false(1, num_triangles);
    
    % Pass the clib objects and the uint64 lengths
    clib.BinaryShader.BinaryRenderer(verticesArg,triangleIDsArg,shadedArg, dir(1),dir(2),dir(3));

    % 5. Convert back to MATLAB to see the result
    ind_gpu_visible(:) = logical(shadedArg);
    if coder.target('MATLAB')
        fprintf('GPU marked %d triangles as visible\n', sum(ind_gpu_visible));
    end
    
    
    ind_visible = ind_rear_facing | ind_gpu_visible;
    ind_shadowed = ~ind_visible;
else
    % All triangles are rear-facing, so all are shadowed
    ind_shadowed = false(1, num_triangles);
end
end

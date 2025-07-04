function [aeroForce__N, aeroTorque__Nm] = newModel(areas__m2,...
    normals,...
    centroids__m,...
    v_rels__m_per_s,...
    deltas__rad,...
    density__kg_per_m3,...
    LUT_data)
    %% newModel - computes aerodynamic forces based on the new IRS Model.
    % Inputs:
    %   areas__m2: 1xN array of the areas of N triangles
    %   normals: 3xN array of surface normals of N triangles
    %   centroids__m: 3xN array of surface centroids of N triangles
    %   v_rels__m_per_s: 3xN array of relative velocities of N triangles
    %   deltas__rad: 1xN array of angles between the flow direction and the normals of N triangles
    %   density__kg_per_m3: scalar value of the incomming streams density
    %   LUT_data: a 'griddedInterpolant' object containing the lookup table data for
    %             the 4 aerodynamic coefficients:
    %             -  C_l_ram
    %             -  C_d_ram
    %             -  C_l_wake
    %             -  C_d_wake
    % Outputs:
    %   aeroForce__N: 3x1 array of the aerodynamic force acting on the body in the same coordinate
    %                 system as the inputs normals and centroids
    %   aeroTorque__Nm: 3x1 array of the aerodynamic torque acting on the body in the same
    %                   coordinate system as the inputs normals and centroids and with respect to its origin 
    %
    arguments
        areas__m2 (1,:) {mustBeNumeric, mustBeReal, mustBePositive};
        normals (3,:) {mustBeNumeric, mustBeReal};
        centroids__m (3,:) {mustBeNumeric, mustBeReal};
        v_rels__m_per_s (3,:) {mustBeNumeric, mustBeReal};
        deltas__rad (1,:) {mustBeNumeric, mustBeReal};
        density__kg_per_m3 (1,1) {mustBeNumeric, mustBeReal, mustBePositive};
        LUT_data {mustBeA(LUT_data, 'griddedInterpolant')};
    end
    %% assert that the return of Lut_data is a nx4 matrix
    assert(isequal(size(LUT_data.Values,2),4), ...
        'LUT_data must return a matrix with 4 columns for C_l_ram, C_d_ram, C_l_wake, C_d_wake');
    %% Abbreviations
    v_rels = v_rels__m_per_s;
    V = vecnorm(v_rels);
    rho = density__kg_per_m3;
    v_hat = v_rels./V;

    %%LUT
    AOA__deg = 90-deltas__rad*180/pi;
    wake_faces = AOA__deg < 0;
    AOA__deg = abs(AOA__deg);
    c = LUT_data(AOA__deg);
    C_l_ram = c(:,1)';
    C_d_ram = c(:,2)';
    C_l_wake = c(:,3)';
    C_d_wake = c(:,4)';

    C_d = C_d_ram;
    C_l = C_l_ram;
    C_d(wake_faces) = C_d_wake(wake_faces);
    C_l(wake_faces) = C_l_wake(wake_faces);


    %darg
    F_d_mag = 0.5*rho*V.^2.*areas__m2.*C_d;
    F_l_mag = 0.5*rho*V.^2.*areas__m2.*C_l;
    F_d = F_d_mag.*v_hat;

    %lift
    lift_dir = -cross(cross(v_hat,normals),v_hat);
    lift_dir_norm = vecnorm(lift_dir);
    lift_dir_norm(lift_dir_norm == 0) = 1; % avoid division by zero
    F_l = F_l_mag.*lift_dir./lift_dir_norm;

    %resultant force
    aeroForce__N = F_l + F_d;
    aeroTorque__Nm = cross(centroids__m,aeroForce__N);
end
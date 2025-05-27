function [aeroForce__N, aeroTorque__Nm] = newModel(areas__m2,...
    normals,...
    centroids__m,...
    v_rels__m_per_s,...
    deltas__rad,...
    density__kg_per_m3,...
    LUT_file)
    %% newModel - computes aerodynamic forces based on the new IRS Model.
    % Inputs:
    %   areas__m2: 1xN array of the areas of N triangles
    %   normals: 3xN array of surface normals of N triangles
    %   centroids__m: 3xN array of surface centroids of N triangles
    %   v_rels__m_per_s: 3xN array of relative velocities of N triangles
    %   deltas__rad: 1xN array of angles between the flow direction and the normals of N triangles
    %   density__kg_per_m3: scalar value of the incomming streams density
    %   LUT_file: path of the .csv-file for the LUT wich holds the c_l and c_d
    %   values depending on the AOA.
    %
    % Outputs:
    %   aeroForce__N: 3x1 array of the aerodynamic force acting on the body in the same coordinate
    %                 system as the inputs normals and centroids
    %   aeroTorque__Nm: 3x1 array of the aerodynamic torque acting on the body in the same
    %                   coordinate system as the inputs normals and centroids and with respect to its origin 
    %
    %% Abbreviations
    v_rels = v_rels__m_per_s;
    V = vecnorm(v_rels);
    rho = density__kg_per_m3;
    v_hat = v_rels./V;

    %%LUT
    AOA__deg = 90-deltas__rad*180/pi;
    backward_facing = AOA__deg < 0;
    AOA__deg(backward_facing) = 0;
    lut = readmatrix(LUT_file);
    AOA_lut = lut(:,1);
    C_l_lut = lut(:,2);
    C_d_lut = lut(:,3);
    C_d = interp1(AOA_lut,C_d_lut,AOA__deg,"linear");
    C_l = interp1(AOA_lut,C_l_lut,AOA__deg,"linear");

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
    F_aero = F_l + F_d;
    F_aero(:,backward_facing) = 0; %backward facing faces have no aerodynamic force
    T_aero = cross(centroids__m,F_aero);
    aeroForce__N = sum(F_aero,2);
    aeroTorque__Nm = sum(T_aero,2);
end





   
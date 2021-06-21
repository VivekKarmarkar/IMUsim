function fh = array_pose_estimate_anatomical_frame
    fh = localfunctions;
end

function [T,R,unit_vectors] = pose_estimate_femoral_frame(current_reconstructed_position_val)
    
    current_side = current_reconstructed_position_val.side;

    r_lat_femoral_epicondyle = current_reconstructed_position_val.(current_side + "_lat_femoral_epicondyle");
    r_med_femoral_epicondyle = current_reconstructed_position_val.(current_side + "_med_femoral_epicondyle");
    r_trochanter = current_reconstructed_position_val.(current_reconstructed_position_val.side + "_trochanter");
    r_origin = 0.5*(r_lat_femoral_epicondyle + r_med_femoral_epicondyle);
    
    ptCloud = pointCloud(transpose(horzcat(r_lat_femoral_epicondyle, r_med_femoral_epicondyle, r_trochanter)));
    pts_collinear = true;
    if ~any(any(isnan(ptCloud.Location)))
        [~,plane_inlierIndices,~] = pcfitplane(ptCloud,1);
        if ~isempty(plane_inlierIndices)
            pts_collinear = false;
        end
    end

    i = [1;0;0];
    j = [0;1;0];
    k = [0;0;1];

    T = r_origin;
    R = nan(3,3);
    R_calc = nan(3,3);
    
    i_hat = nan(3,1);
    j_hat = nan(3,1);
    k_hat = nan(3,1);

    if ~pts_collinear
        
        ry = r_trochanter - r_origin;
        
        r_plane = r_lat_femoral_epicondyle - r_med_femoral_epicondyle;
        r_normal = cross(r_plane, ry);
        
        syms u_x u_y u_z;
        eqn_1 = r_normal(1)*u_x + r_normal(2)*u_y + r_normal(3)*u_z == 0;
        eqn_2 = ry(1)*u_x + ry(2)*u_y + ry(3)*u_z == 0;
        eqn_3 = u_x^2 + u_y^2 + u_z^2 == 1;
        [sol_x, sol_y, sol_z] = solve(eqn_1, eqn_2, eqn_3);
        rz_1 = [double(sol_x(1,1)); double(sol_y(1,1)); double(sol_z(1,1))];
        rz_2 = [double(sol_x(2,1)); double(sol_y(2,1)); double(sol_z(2,1))];
        
        if dot(rz_1, r_plane) < 0
            rz_minus = rz_1;
            rz_plus = rz_2;
        else
            rz_minus = rz_2;
            rz_plus = rz_1;
        end
        
        if current_side == "left"
            rz = rz_minus;
        else
            rz = rz_plus;
        end
        
        rx = cross(ry, rz);
        
        i_hat_dir = rx;
        j_hat_dir = ry;
        k_hat_dir = rz;
        
        i_hat = i_hat_dir/norm(i_hat_dir);
        j_hat = j_hat_dir/norm(j_hat_dir);
        k_hat = k_hat_dir/norm(k_hat_dir);

        R_calc = [dot(i_hat,i) dot(j_hat,i) dot(k_hat,i)
             dot(i_hat,j) dot(j_hat,j) dot(k_hat,j)
             dot(i_hat,k) dot(j_hat,k) dot(k_hat,k)];
    end
    
    tol = 1.0e-10;
    if ~any(any(isnan(R_calc)))
        same_inverse_transpose_bool = all(all(abs(inv(R_calc) - R_calc') < tol ));
        det_R_bool = abs(det(R_calc)-1.0) < tol;
        rotation_matrix_properties_satisfied_bool = same_inverse_transpose_bool && det_R_bool;
        if rotation_matrix_properties_satisfied_bool
            R = R_calc;
        end
    end
    
    unit_vectors.x = i_hat;
    unit_vectors.y = j_hat;
    unit_vectors.z = k_hat;
end

function [T,R,unit_vectors] = pose_estimate_tibial_frame(current_reconstructed_position_val)

    current_side = current_reconstructed_position_val.side;
    r_fibular_head = current_reconstructed_position_val.(current_side + "_fibular_head");
    r_tibial_tuberosity = current_reconstructed_position_val.(current_side + "_tibial_tuberosity");
    r_lateral_malleolus = current_reconstructed_position_val.(current_side + "_lateral_malleolus");
    r_medial_malleolus = current_reconstructed_position_val.(current_side + "_medial_malleolus");
    
    r_origin = 0.5*(r_lateral_malleolus + r_medial_malleolus);
       
    ptCloud = pointCloud(transpose(horzcat(r_fibular_head, r_tibial_tuberosity, r_lateral_malleolus, r_medial_malleolus)));
    pts_collinear = true;
    if ~any(any(isnan(ptCloud.Location)))
        [~,plane_inlierIndices,~] = pcfitplane(ptCloud,1);
        if ~isempty(plane_inlierIndices)
            pts_collinear = false;
        end
    end

    i = [1;0;0];
    j = [0;1;0];
    k = [0;0;1];

    T = r_origin;
    R = nan(3,3);
    R_calc = nan(3,3);
    
    i_hat = nan(3,1);
    j_hat = nan(3,1);
    k_hat = nan(3,1);

    if ~pts_collinear
        
        r_plane1_vector1 = r_lateral_malleolus - r_medial_malleolus;
        r_plane1_vector2 = r_fibular_head - r_origin;
        r_plane1_normal = cross(r_plane1_vector1, r_plane1_vector2);
        
        r_plane2_vector1 = r_plane1_normal;
        r_plane2_vector2 = r_tibial_tuberosity - r_origin;
        r_plane2_normal = cross(r_plane2_vector1, r_plane2_vector2);
        
        ry_1 = cross(r_plane1_normal, r_plane2_normal);
        ry_2 = -ry_1;
        
        if dot(ry_1, r_plane1_vector2) > 0
            ry = ry_1;
        else
            ry = ry_2;
        end
        
        syms u_x u_y u_z;
        eqn_1 = r_plane1_normal(1)*u_x + r_plane1_normal(2)*u_y + r_plane1_normal(3)*u_z == 0;
        eqn_2 = ry(1)*u_x + ry(2)*u_y + ry(3)*u_z == 0;
        eqn_3 = u_x^2 + u_y^2 + u_z^2 == 1;
        [sol_x, sol_y, sol_z] = solve(eqn_1, eqn_2, eqn_3);
        rz_1 = [double(sol_x(1,1)); double(sol_y(1,1)); double(sol_z(1,1))];
        rz_2 = [double(sol_x(2,1)); double(sol_y(2,1)); double(sol_z(2,1))];
        
        if dot(rz_1, r_plane1_vector1) < 0
            rz_minus = rz_1;
            rz_plus = rz_2;
        else
            rz_minus = rz_2;
            rz_plus = rz_1;
        end
        
        if current_side == "left"
            rz = rz_minus;
        else
            rz = rz_plus;
        end
        
        rx = cross(ry, rz);
        
        i_hat_dir = rx;
        j_hat_dir = ry;
        k_hat_dir = rz;
        
        i_hat = i_hat_dir/norm(i_hat_dir);
        j_hat = j_hat_dir/norm(j_hat_dir);
        k_hat = k_hat_dir/norm(k_hat_dir);

        R_calc = [dot(i_hat,i) dot(j_hat,i) dot(k_hat,i)
             dot(i_hat,j) dot(j_hat,j) dot(k_hat,j)
             dot(i_hat,k) dot(j_hat,k) dot(k_hat,k)];
    end
    
    tol = 1.0e-10;
    if ~any(any(isnan(R_calc)))
        same_inverse_transpose_bool = all(all(abs(inv(R_calc) - R_calc') < tol ));
        det_R_bool = abs(det(R_calc)-1.0) < tol;
        rotation_matrix_properties_satisfied_bool = same_inverse_transpose_bool && det_R_bool;
        if rotation_matrix_properties_satisfied_bool
            R = R_calc;
        end
    end
    
    unit_vectors.x = i_hat;
    unit_vectors.y = j_hat;
    unit_vectors.z = k_hat;
end
function fh = array_pose_estimate_anatomical_frame
    fh = localfunctions;
end

function [T,R,unit_vectors] = pose_estimate_femoral_frame(ptCloud)
    pts_collinear = true;

    [~,plane_inlierIndices,~] = pcfitplane(ptCloud,1);
    if ~isempty(plane_inlierIndices)
        pts_collinear = false;
    end

    i = [1;0;0];
    j = [0;1;0];
    k = [0;0;1];
    
    r1 = ptCloud.Location(:,1);
    r2 = ptCloud.Location(:,2);
    r3 = ptCloud.Location(:,3);

    T = mean(horzcat(r1,r2,r3),2);
    R = nan(3,3);
    R_calc = nan(3,3);
    
    i_hat = nan(3,1);
    j_hat = nan(3,1);
    k_hat = nan(3,1);

    if ~pts_collinear

        r21 = r2 - r1;
        r31 = r3 - r1;

        u21 = r21/norm(r21);
        u31 = r31/norm(r31);

        i_hat_dir = u21;
        j_hat_dir = cross(u21, u31);
        k_hat_dir = cross(i_hat_dir, j_hat_dir);
        
        i_hat = i_hat_dir/norm(i_hat_dir);
        j_hat = j_hat_dir/norm(j_hat_dir);
        k_hat = k_hat_dir/norm(k_hat_dir);

        R_calc = [dot(i_hat,i) dot(j_hat,i) dot(k_hat,i)
             dot(i_hat,j) dot(j_hat,j) dot(k_hat,j)
             dot(i_hat,k) dot(j_hat,k) dot(k_hat,k)];
    end
    
    tol = 1.0e-10;
    same_inverse_transpose_bool = all(all(abs(inv(R_calc) - R_calc') < tol ));
    det_R_bool = abs(det(R_calc)-1.0) < tol;
    rotation_matrix_properties_satisfied_bool = same_inverse_transpose_bool && det_R_bool;
    if rotation_matrix_properties_satisfied_bool
        R = R_calc;
    end
    
    unit_vectors.x = i_hat;
    unit_vectors.y = j_hat;
    unit_vectors.z = k_hat;
end

function [T,R,unit_vectors] = pose_estimate_tibial_frame(ptCloud)
    pts_collinear = true;

    [~,plane_inlierIndices,~] = pcfitplane(ptCloud,1);
    if ~isempty(plane_inlierIndices)
        pts_collinear = false;
    end

    i = [1;0;0];
    j = [0;1;0];
    k = [0;0;1];
    
    r1 = ptCloud.Location(:,1);
    r2 = ptCloud.Location(:,2);
    r3 = ptCloud.Location(:,3);

    T = mean(horzcat(r1,r2,r3),2);
    R = nan(3,3);
    R_calc = nan(3,3);
    
    i_hat = nan(3,1);
    j_hat = nan(3,1);
    k_hat = nan(3,1);

    if ~pts_collinear

        r21 = r2 - r1;
        r31 = r3 - r1;

        u21 = r21/norm(r21);
        u31 = r31/norm(r31);

        i_hat_dir = u21;
        j_hat_dir = cross(u21, u31);
        k_hat_dir = cross(i_hat_dir, j_hat_dir);
        
        i_hat = i_hat_dir/norm(i_hat_dir);
        j_hat = j_hat_dir/norm(j_hat_dir);
        k_hat = k_hat_dir/norm(k_hat_dir);

        R_calc = [dot(i_hat,i) dot(j_hat,i) dot(k_hat,i)
             dot(i_hat,j) dot(j_hat,j) dot(k_hat,j)
             dot(i_hat,k) dot(j_hat,k) dot(k_hat,k)];
    end
    
    tol = 1.0e-10;
    same_inverse_transpose_bool = all(all(abs(inv(R_calc) - R_calc') < tol ));
    det_R_bool = abs(det(R_calc)-1.0) < tol;
    rotation_matrix_properties_satisfied_bool = same_inverse_transpose_bool && det_R_bool;
    if rotation_matrix_properties_satisfied_bool
        R = R_calc;
    end
    
    unit_vectors.x = i_hat;
    unit_vectors.y = j_hat;
    unit_vectors.z = k_hat;
end
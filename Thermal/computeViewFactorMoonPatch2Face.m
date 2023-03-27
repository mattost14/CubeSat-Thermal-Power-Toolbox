function F = computeViewFactorMoonPatch2Face(r_sat, v_sat, i_patch, R, param)

    % function that return angle between two vectors in degrees
    calculateAngleBetweenVectors = @(u,v) acosd(max(min(dot(u,v)/(norm(u)*norm(v)),1),-1));

    % Get face vector at LVLH frame
    XplusFaceVector_LVLH = param.pwr.XplusFaceVector_LVLH;
    YplusFaceVector_LVLH = param.pwr.YplusFaceVector_LVLH;
    ZplusFaceVector_LVLH = param.pwr.ZplusFaceVector_LVLH;

    % Get rotation matrix from LVLH to Fixed Frame
    z_L2F = -r_sat/norm(r_sat);   % Unit vector z
    h_F = cross(r_sat, v_sat);    % Angular momentum 
    y_L2F = -h_F/norm(h_F);             % Unit vector y
    x_L2F = cross(y_L2F, z_L2F);        % Unit vector x
    A_F2L = [x_L2F'; y_L2F'; z_L2F'];   % Rotation matrix from Fixed Frame to LVLH frame
    A_L2F = A_F2L';                     % Rotation matrix from LVLH frame to Fixed Frame
    
    % Get face vector at Fixed Frame
    XplusFaceVector_F = A_L2F * XplusFaceVector_LVLH;
    YplusFaceVector_F = A_L2F * YplusFaceVector_LVLH;
    ZplusFaceVector_F = A_L2F * ZplusFaceVector_LVLH;

    % View factor function (ref:
    % http://imartinez.etsiae.upm.es/~isidoro/tc3/Radiation%20View%20factors.pdf,
    % Eq. 1, Pag. 3)
    computeViewFactor = @(beta1, beta2, A2, r12) cosd(beta1)*cosd(beta2)*A2/(pi*r12^2);

    % ### Compute View Factor ###
    patch2sat = r_sat - i_patch * R;
    i_patch2sat = patch2sat/norm(patch2sat);
    beta1 = calculateAngleBetweenVectors(i_patch, patch2sat);

    beta2(1) =  calculateAngleBetweenVectors(XplusFaceVector_F, -i_patch2sat);
    beta2(2) =  calculateAngleBetweenVectors(-XplusFaceVector_F, -i_patch2sat);
    beta2(3) =  calculateAngleBetweenVectors(YplusFaceVector_F, -i_patch2sat);
    beta2(4) =  calculateAngleBetweenVectors(-YplusFaceVector_F, -i_patch2sat);
    beta2(5) =  calculateAngleBetweenVectors(ZplusFaceVector_F, -i_patch2sat);
    beta2(6) =  calculateAngleBetweenVectors(-ZplusFaceVector_F, -i_patch2sat);

    F = zeros(6,1);
    for n=1:6
        if(beta2(n)<90)
            F(n) = computeViewFactor(beta1, beta2(n), param.facesArea(n), norm(patch2sat));
        end
    end
end

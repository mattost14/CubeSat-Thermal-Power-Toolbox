function plotCubesat(param)
    facesPoints = getCubeSatPlotPoints(param.numberOfUnits);
    X = [1 0 0];
    Y = [0 1 0];
    Z = [0 0 1];
    scale = 2.9;
    for i=1:length(param.facesMaterial)
        h = surf(facesPoints(i).x, facesPoints(i).y, facesPoints(i).z, 'CData',getTexture(param.facesMaterial(i)), 'FaceColor','texturemap', 'edgecolor','#656565', 'LineWidth', 3); 
%         rotate(h, Z, param.attitude.rotationAboutZ, [0,0,0]) %rotation about Z
%         rotate(h, angle2dcm(deg2rad(param.attitude.rotationAboutZ), 0, 0,'ZYX')'*Y', param.attitude.rotationAboutY, [0,0,0]) %rotation about Y
%         rotate(h, angle2dcm(deg2rad(param.attitude.rotationAboutZ), deg2rad(param.attitude.rotationAboutY), 0,'ZYX')'*X', param.attitude.rotationAboutX, [0,0,0]) %rotation about X
        hold on;
    end
    
    plot3([-5,5],[0,0],[0,0],'LineWidth',2, 'Color', '#f2b009')
    
    h = surf([-3, 3; -3, 3],[-3, -3; 3, 3], [-3, -3; -3, -3], 'CData',getTexture(param.centralBody), 'FaceColor','texturemap', 'edgecolor','#656565', 'LineWidth', 3);

    axis([ -1  1  -1  1 -1  1]*3)
    grid on
    axis vis3d
    
end
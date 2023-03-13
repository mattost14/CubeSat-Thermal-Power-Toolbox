function  plotOrbit3D(appAxis, centralBody, r_ff)
%PLOT3DORBIT Summary of this function goes here
%   Detailed explanation goes here
    cla(appAxis);
    switch centralBody
        case "Earth"
            C = imread('Textures/earthGlobe.png');
            R = 6378.135;
            [x, y, z] = ellipsoid(0,0,0, R, R, R,1E2);
            surf(appAxis, x,y,z,circshift(flip(C),[0,ceil(size(C,2)/360)]), 'FaceColor', 'texturemap','EdgeColor','none');
        case "Moon"
            C = imread('Textures/moonGlobe.png','png');
            if(size(C,3)==1)
                C = cat(3, C, C, C); 
            end
            R = 1737.4;
            [x, y, z] = ellipsoid(0,0,0, R, R, R,1E2);
            surf(appAxis, x,y,z,circshift(flip(C),[0,ceil(size(C,2)/360)]), 'FaceColor', 'texturemap','EdgeColor','none');
    end

    hold(appAxis, 'on')
    v = 2*R*[1;0;0];
    q=quiver3(appAxis, [0],[0],[0],v(1),v(2),v(3), 'Color','g','LineWidth',5);
    text(appAxis, v(1),v(2),v(3),'x','HorizontalAlignment','left','FontSize',12, 'Color','g');
    v = 2*R*[0;1;0];
    q=quiver3(appAxis, [0],[0],[0],v(1),v(2),v(3), 'Color','g','LineWidth',5);
    text(appAxis, v(1),v(2),v(3),'y','HorizontalAlignment','left','FontSize',12, 'Color','r');
    v = 2*R*[0;0;1];
    q=quiver3(appAxis, [0],[0],[0],v(1),v(2),v(3), 'Color','g','LineWidth',5);
    text(appAxis, v(1),v(2),v(3),'z','HorizontalAlignment','left','FontSize',12, 'Color','b');

    
    axis(appAxis, "off");
    axis(appAxis, "equal");

    if(nargin>2)
        plot3(appAxis, r_ff(:,1)/1e3, r_ff(:,2)/1e3, r_ff(:,3)/1e3,"r");
    end
end


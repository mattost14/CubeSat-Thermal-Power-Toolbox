function  plotOrbit3D(appAxis, centralBody, param)
%PLOT3DORBIT plots the propagated orbit around the selected central body in
%the Fixed Frame

    
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
        % Check if there orbital was propagated successfully
        if(param.orb.prop.flag)
            r_ff = param.orb.prop.r_ff;
            plot3(appAxis, r_ff(:,1)/1e3, r_ff(:,2)/1e3, r_ff(:,3)/1e3,"r");
            scatter3(appAxis, r_ff(1,1)/1e3, r_ff(1,2)/1e3, r_ff(1,3)/1e3,"r",'filled'); % Show sat dot at initial position
        end
        % Check if the user wants to show the Sun vector at time t=0
        if(param.orb.prop.flag && param.plot.orbitPlot.showSunVectorFlag)
            sunVector = param.orb.prop.Sun_ff;
            v = 2*R*sunVector(1,:); % First sun position in the Fixed Frame
            quiver3(appAxis, [0],[0],[0],v(1),v(2),v(3), 'Color',"#EDB120",'LineWidth',5);
            text(appAxis, v(1),v(2),v(3),'Sun','HorizontalAlignment','left','FontSize',12, 'Color',"#EDB120");  
        end
    end
end


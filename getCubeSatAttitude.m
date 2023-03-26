function [XplusFaceVector_LVLH, YplusFaceVector_LVLH, ZplusFaceVector_LVLH] = getCubeSatAttitude(param)
    % Get attitude settings
    nadirFace = param.attitude.nadirFace;
    ramFace = param.attitude.ramFace;
    
    % Define sat attitude in the LVLH based on the attitude settings
    nadirFaceId=char(nadirFace); 
    switch nadirFaceId(1) 
        case "X"
            if(nadirFaceId(2)=="+")
                XplusFaceVector_LVLH = [0; 0; 1];
            else
                XplusFaceVector_LVLH = [0; 0; -1];
            end
        case "Y"
            if(nadirFaceId(2)=="+")
                YplusFaceVector_LVLH = [0; 0; 1];
            else
                YplusFaceVector_LVLH = [0; 0; -1];
            end
        case "Z"
            if(nadirFaceId(2)=="+")
                ZplusFaceVector_LVLH = [0; 0; 1];
            else
                ZplusFaceVector_LVLH = [0; 0; -1];
            end
    end
    
    ramFaceId=char(ramFace); 
    switch ramFaceId(1) 
        case "X"
            if(ramFaceId(2)=="+")
                XplusFaceVector_LVLH = [1; 0; 0];
            else
                XplusFaceVector_LVLH = [-1; 0; 0];
            end
        case "Y"
            if(ramFaceId(2)=="+")
                YplusFaceVector_LVLH = [1; 0; 0];
            else
                YplusFaceVector_LVLH = [-1; 0; 0];
            end
        case "Z"
            if(ramFaceId(2)=="+")
                ZplusFaceVector_LVLH = [1; 0; 0];
            else
                ZplusFaceVector_LVLH = [-1; 0; 0];
            end
    end
    
    listOfFaces = ["X", "Y", "Z"];
    filterCondition = find((listOfFaces~=nadirFaceId(1)) & (listOfFaces~=ramFaceId(1)));
    lastFace = listOfFaces(filterCondition);
    switch lastFace
        case "X"
                XplusFaceVector_LVLH = cross(YplusFaceVector_LVLH, ZplusFaceVector_LVLH);
        case "Y"
                YplusFaceVector_LVLH = cross(ZplusFaceVector_LVLH, XplusFaceVector_LVLH);
        case "Z"
                ZplusFaceVector_LVLH = cross(XplusFaceVector_LVLH, YplusFaceVector_LVLH);
    end
end
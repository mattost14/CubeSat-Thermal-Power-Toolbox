function [facesBlanketMass, resistanceBlanket, resistanceInterFacesBlanket] = setBlanketParameters(param)
%SETBLANKETPARAMETERS sets up the thermal parameters in case the face is covered with
% MLI blanket

    facesBlanketMass = zeros(1,6);
    resistanceBlanket = zeros(1,6);
    resistanceInterFacesBlanket = ones(6,6)*1e3;

    blanketOptions = getPropertiesFromMaterial("blankets");
    faceText = ["X+","X-","Y+","Y-","Z+","Z-"];
    for n=1:6
        if(any(strcmp(blanketOptions,param.facesMaterial(n)))) % if the face material is a "blanket" (insulator)
            face = char(faceText(n));
            % Get blanket properties
            p = getPropertiesFromMaterial(param.facesMaterial(n));
            % Adjust the face blanket mass
            facesBlanketMass(n) = p.areaDensity * param.facesArea(n); % kg
            % Adjust the blanket resistance
            resistanceBlanket(n) = p.conductionResistanceTimesArea / param.facesArea(n);    % K/W (High thermal resistance crossing the blanket)
            
            % Find neighboor faces that is also covered with the MLI
            % blanket
            for otherFacesID=1:6
                otherFace = char(faceText(otherFacesID));
                % If the neighboor face is also covered with MLI - Blanket
                if(otherFace(1)~=face(1) && any(strcmp(blanketOptions, param.facesMaterial(otherFacesID))) )
                    resistanceInterFacesBlanket(n,otherFacesID) = .1;        % K/W (Low thermal resistance along the blanket surface)
                elseif(otherFace(1)==face(1))
                    resistanceInterFacesBlanket(n,otherFacesID) = nan;        % Same or oposite face (Not applicable)
                else
                    resistanceInterFacesBlanket(n,otherFacesID) = 1e6;       % K/W (High thermal resistance)
                end
            end
        else
            resistanceBlanket(n) = 0;    % K/W (No blanket -> no resistance)
            facesBlanketMass(n) = 1e-5;     % kg (No mass)
        end
    end
end


function thermalModel = prepareSimscapeModelInputs_v2(param)
    % ######## PREPARE SIMSCAPE simpleThermalModel MODEL ########
    thermalModel.mdl = "simpleThermalModel";

    % Use a SimulationInput object to configure the model for our simulation.
    thermalModel.sim.in = Simulink.SimulationInput(thermalModel.mdl);

    % Solver setttings
    thermalModel.sim.in = thermalModel.sim.in.setModelParameter(...
        SolverType="Variable-step",...
        SolverName="VariableStepAuto",...
        RelTol="1e-6",...
        AbsTol="1e-7",...
        StopTime=string(param.orb.durationHours*3600));
    
    
    % Define thermal resistance inputs
    blk_RXplus = thermalModel.mdl + "/RX+";
    blk_RXminus = thermalModel.mdl + "/RX-";
    blk_RYplus = thermalModel.mdl + "/RY+";
    blk_RYminus = thermalModel.mdl + "/RY-";
    blk_RZplus = thermalModel.mdl + "/RZ+";
    blk_RZminus = thermalModel.mdl + "/RZ-";

    thermalModel.sim.in = thermalModel.sim.in.setBlockParameter(...
    blk_RXplus, "resistance", string(param.resistance(1)),...
    blk_RXminus, "resistance", string(param.resistance(2)),...
    blk_RYplus, "resistance", string(param.resistance(3)),...
    blk_RYminus, "resistance", string(param.resistance(4)),...
    blk_RZplus, "resistance", string(param.resistance(5)),...
    blk_RZminus, "resistance", string(param.resistance(6)));

    %% Define inputs for the Internal Node 1
    blk_thermalMass = thermalModel.mdl + "/Internal Node/Internal Node 1";
    blk_internalNodePwr = thermalModel.mdl + "/Internal Node/Internal Node Pwr";

    thermalModel.sim.in = thermalModel.sim.in.setBlockParameter(...
    blk_thermalMass, "mass", string(param.internalNodeMass(1)),...
    blk_thermalMass, "sp_heat", string(param.internalNodeCp(1)),...
    blk_internalNodePwr, "VariableName", 'param.thermal.internalPwr');


    %% Define inputs for Face X+ block
    faceId = 1;

    % Main
    blk_RadBlanket = thermalModel.mdl + "/Face X+/RadHeatBlanket";
    blk_Rblanket = thermalModel.mdl + "/Face X+/R_Blanket";
    blk_thermalMass = thermalModel.mdl + "/Face X+/Thermal Mass";

    thermalModel.sim.in = thermalModel.sim.in.setBlockParameter(... 
        blk_RadBlanket, "area", string(param.facesArea(faceId)),...
        blk_RadBlanket, "rad_tr_coeff", string(param.kradBlanket(faceId)),...
        blk_Rblanket, "resistance", string(param.resistanceBlanket(faceId)),...
        blk_thermalMass, "mass", string(param.facesMass(faceId)),...
        blk_thermalMass, "sp_heat", string(param.facesCp(faceId)));

    % Enviromental Block
    blk_solarFlux =   thermalModel.mdl + "/Thermal Loads from Enviroment X+/Solar Flux";
    blk_IRFlux =   thermalModel.mdl + "/Thermal Loads from Enviroment X+/IR Flux";
    blk_AlbedoFlux =   thermalModel.mdl + "/Thermal Loads from Enviroment X+/Albedo Flux";
    
    thermalModel.sim.in = thermalModel.sim.in.setBlockParameter(...   
        blk_solarFlux, "VariableName", 'param.thermal.faceSolar.Xplus',...
        blk_IRFlux, "VariableName", 'param.thermal.faceIR.Xplus',...
        blk_AlbedoFlux, "VariableName",'param.thermal.faceAlbedo.Xplus');

    % Neighboor Faces Sub-block
    blk_RXplusYplus = thermalModel.mdl + "/Face X+/Neighboor Faces/RX+Y+";
    blk_RXplusYminus = thermalModel.mdl + "/Face X+/Neighboor Faces/RX+Y-";
    blk_RXplusZplus = thermalModel.mdl + "/Face X+/Neighboor Faces/RX+Z+";
    blk_RXplusZminus = thermalModel.mdl + "/Face X+/Neighboor Faces/RX+Z-";

    thermalModel.sim.in = thermalModel.sim.in.setBlockParameter(...
        blk_RXplusYplus, "resistance", string(param.resistanceInterFaces(faceId, 3)),...
        blk_RXplusYminus, "resistance", string(param.resistanceInterFaces(faceId, 4)),...
        blk_RXplusZplus, "resistance", string(param.resistanceInterFaces(faceId, 5)),...
        blk_RXplusZminus, "resistance", string(param.resistanceInterFaces(faceId, 6)));

    % OuterLayer Sub-block
    blk_radiativeHeatTransfer = thermalModel.mdl + "/Face X+/OuterLayer/Radiative Heat Transfer";
    blk_thermalMassBlanket = thermalModel.mdl + "/Face X+/OuterLayer/Thermal Mass Blanket";


    thermalModel.sim.in = thermalModel.sim.in.setBlockParameter(...   
        blk_radiativeHeatTransfer, "area", string(param.facesArea(faceId)),...
        blk_radiativeHeatTransfer, "rad_tr_coeff", string(param.krad(faceId)),...
        blk_thermalMassBlanket, "mass", string(param.facesBlanketMass(faceId)),...
        blk_thermalMassBlanket, "sp_heat", string(param.facesBlanketCp(faceId)));

    %% Define inputs for Face X- block
    faceId = 2;

    % Main
    blk_RadBlanket = thermalModel.mdl + "/Face X-/RadHeatBlanket";
    blk_Rblanket = thermalModel.mdl + "/Face X-/R_Blanket";
    blk_thermalMass = thermalModel.mdl + "/Face X-/Thermal Mass";

    thermalModel.sim.in = thermalModel.sim.in.setBlockParameter(... 
        blk_RadBlanket, "area", string(param.facesArea(faceId)),...
        blk_RadBlanket, "rad_tr_coeff", string(param.kradBlanket(faceId)),...
        blk_Rblanket, "resistance", string(param.resistanceBlanket(faceId)),...
        blk_thermalMass, "mass", string(param.facesMass(faceId)),...
        blk_thermalMass, "sp_heat", string(param.facesCp(faceId)));

    % Enviromental Block
    blk_solarFlux =   thermalModel.mdl + "/Thermal Loads from Enviroment X-/Solar Flux";
    blk_IRFlux =   thermalModel.mdl + "/Thermal Loads from Enviroment X-/IR Flux";
    blk_AlbedoFlux =   thermalModel.mdl + "/Thermal Loads from Enviroment X-/Albedo Flux";
    
    thermalModel.sim.in = thermalModel.sim.in.setBlockParameter(...   
        blk_solarFlux, "VariableName", 'param.thermal.faceSolar.Xminus',...
        blk_IRFlux, "VariableName", 'param.thermal.faceIR.Xminus',...
        blk_AlbedoFlux, "VariableName",'param.thermal.faceAlbedo.Xminus');

    % Neighboor Faces Sub-block
    blk_RXminusYplus = thermalModel.mdl + "/Face X-/Neighboor Faces/RX-Y+";
    blk_RXminusYminus = thermalModel.mdl + "/Face X-/Neighboor Faces/RX-Y-";
    blk_RXminusZplus = thermalModel.mdl + "/Face X-/Neighboor Faces/RX-Z+";
    blk_RXminusZminus = thermalModel.mdl + "/Face X-/Neighboor Faces/RX-Z-";

    thermalModel.sim.in = thermalModel.sim.in.setBlockParameter(...
        blk_RXminusYplus, "resistance", string(param.resistanceInterFaces(faceId, 3)),...
        blk_RXminusYminus, "resistance", string(param.resistanceInterFaces(faceId, 4)),...
        blk_RXminusZplus, "resistance", string(param.resistanceInterFaces(faceId, 5)),...
        blk_RXminusZminus, "resistance", string(param.resistanceInterFaces(faceId, 6)));

    % OuterLayer Sub-block
    blk_radiativeHeatTransfer = thermalModel.mdl + "/Face X-/OuterLayer/Radiative Heat Transfer";
    blk_thermalMassBlanket = thermalModel.mdl + "/Face X-/OuterLayer/Thermal Mass Blanket";


    thermalModel.sim.in = thermalModel.sim.in.setBlockParameter(...   
        blk_radiativeHeatTransfer, "area", string(param.facesArea(faceId)),...
        blk_radiativeHeatTransfer, "rad_tr_coeff", string(param.krad(faceId)),...
        blk_thermalMassBlanket, "mass", string(param.facesBlanketMass(faceId)),...
        blk_thermalMassBlanket, "sp_heat", string(param.facesBlanketCp(faceId)));

    %% Define inputs for Face Y+ block
    faceId = 3;

    % Main
    blk_RadBlanket = thermalModel.mdl + "/Face Y+/RadHeatBlanket";
    blk_Rblanket = thermalModel.mdl + "/Face Y+/R_Blanket";
    blk_thermalMass = thermalModel.mdl + "/Face Y+/Thermal Mass";

    thermalModel.sim.in = thermalModel.sim.in.setBlockParameter(... 
        blk_RadBlanket, "area", string(param.facesArea(faceId)),...
        blk_RadBlanket, "rad_tr_coeff", string(param.kradBlanket(faceId)),...
        blk_Rblanket, "resistance", string(param.resistanceBlanket(faceId)),...
        blk_thermalMass, "mass", string(param.facesMass(faceId)),...
        blk_thermalMass, "sp_heat", string(param.facesCp(faceId)));

    % Enviromental Block
    blk_solarFlux =   thermalModel.mdl + "/Thermal Loads from Enviroment Y+/Solar Flux";
    blk_IRFlux =   thermalModel.mdl + "/Thermal Loads from Enviroment Y+/IR Flux";
    blk_AlbedoFlux =   thermalModel.mdl + "/Thermal Loads from Enviroment Y+/Albedo Flux";
    
    thermalModel.sim.in = thermalModel.sim.in.setBlockParameter(...   
        blk_solarFlux, "VariableName", 'param.thermal.faceSolar.Yplus',...
        blk_IRFlux, "VariableName", 'param.thermal.faceIR.Yplus',...
        blk_AlbedoFlux, "VariableName",'param.thermal.faceAlbedo.Yplus');

    % Neighboor Faces Sub-block
    blk_RYplusXplus = thermalModel.mdl + "/Face Y+/Neighboor Faces/RY+X+";
    blk_RYplusXminus = thermalModel.mdl + "/Face Y+/Neighboor Faces/RY+X-";
    blk_RYplusZplus = thermalModel.mdl + "/Face Y+/Neighboor Faces/RY+Z+";
    blk_RYplusZminus = thermalModel.mdl + "/Face Y+/Neighboor Faces/RY+Z-";

    thermalModel.sim.in = thermalModel.sim.in.setBlockParameter(...
        blk_RYplusXplus, "resistance", string(param.resistanceInterFaces(faceId, 1)),...
        blk_RYplusXminus, "resistance", string(param.resistanceInterFaces(faceId, 2)),...
        blk_RYplusZplus, "resistance", string(param.resistanceInterFaces(faceId, 5)),...
        blk_RYplusZminus, "resistance", string(param.resistanceInterFaces(faceId, 6)));

    % OuterLayer Sub-block
    blk_radiativeHeatTransfer = thermalModel.mdl + "/Face Y+/OuterLayer/Radiative Heat Transfer";
    blk_thermalMassBlanket = thermalModel.mdl + "/Face Y+/OuterLayer/Thermal Mass Blanket";


    thermalModel.sim.in = thermalModel.sim.in.setBlockParameter(...   
        blk_radiativeHeatTransfer, "area", string(param.facesArea(faceId)),...
        blk_radiativeHeatTransfer, "rad_tr_coeff", string(param.krad(faceId)),...
        blk_thermalMassBlanket, "mass", string(param.facesBlanketMass(faceId)),...
        blk_thermalMassBlanket, "sp_heat", string(param.facesBlanketCp(faceId)));

    %% Define inputs for Face Y- block
    faceId = 4;

    % Main
    blk_RadBlanket = thermalModel.mdl + "/Face Y-/RadHeatBlanket";
    blk_Rblanket = thermalModel.mdl + "/Face Y-/R_Blanket";
    blk_thermalMass = thermalModel.mdl + "/Face Y-/Thermal Mass";

    thermalModel.sim.in = thermalModel.sim.in.setBlockParameter(... 
        blk_RadBlanket, "area", string(param.facesArea(faceId)),...
        blk_RadBlanket, "rad_tr_coeff", string(param.kradBlanket(faceId)),...
        blk_Rblanket, "resistance", string(param.resistanceBlanket(faceId)),...
        blk_thermalMass, "mass", string(param.facesMass(faceId)),...
        blk_thermalMass, "sp_heat", string(param.facesCp(faceId)));

    % Enviromental Block
    blk_solarFlux =   thermalModel.mdl + "/Thermal Loads from Enviroment Y-/Solar Flux";
    blk_IRFlux =   thermalModel.mdl + "/Thermal Loads from Enviroment Y-/IR Flux";
    blk_AlbedoFlux =   thermalModel.mdl + "/Thermal Loads from Enviroment Y-/Albedo Flux";
    
    thermalModel.sim.in = thermalModel.sim.in.setBlockParameter(...   
        blk_solarFlux, "VariableName", 'param.thermal.faceSolar.Yminus',...
        blk_IRFlux, "VariableName", 'param.thermal.faceIR.Yminus',...
        blk_AlbedoFlux, "VariableName",'param.thermal.faceAlbedo.Yminus');

    % Neighboor Faces Sub-block
    blk_RYminusXplus = thermalModel.mdl + "/Face Y-/Neighboor Faces/RY-X+";
    blk_RYminusXminus = thermalModel.mdl + "/Face Y-/Neighboor Faces/RY-X-";
    blk_RYminusZplus = thermalModel.mdl + "/Face Y-/Neighboor Faces/RY-Z+";
    blk_RYminusZminus = thermalModel.mdl + "/Face Y-/Neighboor Faces/RY-Z-";

    thermalModel.sim.in = thermalModel.sim.in.setBlockParameter(...
        blk_RYminusXplus, "resistance", string(param.resistanceInterFaces(faceId, 1)),...
        blk_RYminusXminus, "resistance", string(param.resistanceInterFaces(faceId, 2)),...
        blk_RYminusZplus, "resistance", string(param.resistanceInterFaces(faceId, 5)),...
        blk_RYminusZminus, "resistance", string(param.resistanceInterFaces(faceId, 6)));

    % OuterLayer Sub-block
    blk_radiativeHeatTransfer = thermalModel.mdl + "/Face Y-/OuterLayer/Radiative Heat Transfer";
    blk_thermalMassBlanket = thermalModel.mdl + "/Face Y-/OuterLayer/Thermal Mass Blanket";


    thermalModel.sim.in = thermalModel.sim.in.setBlockParameter(...   
        blk_radiativeHeatTransfer, "area", string(param.facesArea(faceId)),...
        blk_radiativeHeatTransfer, "rad_tr_coeff", string(param.krad(faceId)),...
        blk_thermalMassBlanket, "mass", string(param.facesBlanketMass(faceId)),...
        blk_thermalMassBlanket, "sp_heat", string(param.facesBlanketCp(faceId)));

    %% Define inputs for Face Z+ block
    faceId = 5;

    % Main
    blk_RadBlanket = thermalModel.mdl + "/Face Z+/RadHeatBlanket";
    blk_Rblanket = thermalModel.mdl + "/Face Z+/R_Blanket";
    blk_thermalMass = thermalModel.mdl + "/Face Z+/Thermal Mass";

    thermalModel.sim.in = thermalModel.sim.in.setBlockParameter(... 
        blk_RadBlanket, "area", string(param.facesArea(faceId)),...
        blk_RadBlanket, "rad_tr_coeff", string(param.kradBlanket(faceId)),...
        blk_Rblanket, "resistance", string(param.resistanceBlanket(faceId)),...
        blk_thermalMass, "mass", string(param.facesMass(faceId)),...
        blk_thermalMass, "sp_heat", string(param.facesCp(faceId)));

    % Enviromental Block
    blk_solarFlux =   thermalModel.mdl + "/Thermal Loads from Enviroment Z+/Solar Flux";
    blk_IRFlux =   thermalModel.mdl + "/Thermal Loads from Enviroment Z+/IR Flux";
    blk_AlbedoFlux =   thermalModel.mdl + "/Thermal Loads from Enviroment Z+/Albedo Flux";
    
    thermalModel.sim.in = thermalModel.sim.in.setBlockParameter(...   
        blk_solarFlux, "VariableName", 'param.thermal.faceSolar.Zplus',...
        blk_IRFlux, "VariableName", 'param.thermal.faceIR.Zplus',...
        blk_AlbedoFlux, "VariableName",'param.thermal.faceAlbedo.Zplus');

    % Neighboor Faces Sub-block
    blk_RZplusXplus = thermalModel.mdl + "/Face Z+/Neighboor Faces/RZ+X+";
    blk_RZplusXminus = thermalModel.mdl + "/Face Z+/Neighboor Faces/RZ+X-";
    blk_RZplusYplus = thermalModel.mdl + "/Face Z+/Neighboor Faces/RZ+Y+";
    blk_RZplusYminus = thermalModel.mdl + "/Face Z+/Neighboor Faces/RZ+Y-";

    thermalModel.sim.in = thermalModel.sim.in.setBlockParameter(...
        blk_RZplusXplus, "resistance", string(param.resistanceInterFaces(faceId, 1)),...
        blk_RZplusXminus, "resistance", string(param.resistanceInterFaces(faceId, 2)),...
        blk_RZplusYplus, "resistance", string(param.resistanceInterFaces(faceId, 3)),...
        blk_RZplusYminus, "resistance", string(param.resistanceInterFaces(faceId, 4)));

    % OuterLayer Sub-block
    blk_radiativeHeatTransfer = thermalModel.mdl + "/Face Z+/OuterLayer/Radiative Heat Transfer";
    blk_thermalMassBlanket = thermalModel.mdl + "/Face Z+/OuterLayer/Thermal Mass Blanket";


    thermalModel.sim.in = thermalModel.sim.in.setBlockParameter(...   
        blk_radiativeHeatTransfer, "area", string(param.facesArea(faceId)),...
        blk_radiativeHeatTransfer, "rad_tr_coeff", string(param.krad(faceId)),...
        blk_thermalMassBlanket, "mass", string(param.facesBlanketMass(faceId)),...
        blk_thermalMassBlanket, "sp_heat", string(param.facesBlanketCp(faceId)));

    %% Define inputs for Face Z- block
    faceId = 6;

    % Main
    blk_RadBlanket = thermalModel.mdl + "/Face Z-/RadHeatBlanket";
    blk_Rblanket = thermalModel.mdl + "/Face Z-/R_Blanket";
    blk_thermalMass = thermalModel.mdl + "/Face Z-/Thermal Mass";

    thermalModel.sim.in = thermalModel.sim.in.setBlockParameter(... 
        blk_RadBlanket, "area", string(param.facesArea(faceId)),...
        blk_RadBlanket, "rad_tr_coeff", string(param.kradBlanket(faceId)),...
        blk_Rblanket, "resistance", string(param.resistanceBlanket(faceId)),...
        blk_thermalMass, "mass", string(param.facesMass(faceId)),...
        blk_thermalMass, "sp_heat", string(param.facesCp(faceId)));

    % Enviromental Block
    blk_solarFlux =   thermalModel.mdl + "/Thermal Loads from Enviroment Z-/Solar Flux";
    blk_IRFlux =   thermalModel.mdl + "/Thermal Loads from Enviroment Z-/IR Flux";
    blk_AlbedoFlux =   thermalModel.mdl + "/Thermal Loads from Enviroment Z-/Albedo Flux";
    
    thermalModel.sim.in = thermalModel.sim.in.setBlockParameter(...   
        blk_solarFlux, "VariableName", 'param.thermal.faceSolar.Yminus',...
        blk_IRFlux, "VariableName", 'param.thermal.faceIR.Yminus',...
        blk_AlbedoFlux, "VariableName",'param.thermal.faceAlbedo.Yminus');

    % Neighboor Faces Sub-block
    blk_RZminusXplus = thermalModel.mdl + "/Face Z-/Neighboor Faces/RZ-X+";
    blk_RZminusXminus = thermalModel.mdl + "/Face Z-/Neighboor Faces/RZ-X-";
    blk_RZminusYplus = thermalModel.mdl + "/Face Z-/Neighboor Faces/RZ-Y+";
    blk_RZminusYminus = thermalModel.mdl + "/Face Z-/Neighboor Faces/RZ-Y-";

    thermalModel.sim.in = thermalModel.sim.in.setBlockParameter(...
        blk_RZminusXplus, "resistance", string(param.resistanceInterFaces(faceId, 1)),...
        blk_RZminusXminus, "resistance", string(param.resistanceInterFaces(faceId, 2)),...
        blk_RZminusYplus, "resistance", string(param.resistanceInterFaces(faceId, 3)),...
        blk_RZminusYminus, "resistance", string(param.resistanceInterFaces(faceId, 4)));

    % OuterLayer Sub-block
    blk_radiativeHeatTransfer = thermalModel.mdl + "/Face Z-/OuterLayer/Radiative Heat Transfer";
    blk_thermalMassBlanket = thermalModel.mdl + "/Face Z-/OuterLayer/Thermal Mass Blanket";


    thermalModel.sim.in = thermalModel.sim.in.setBlockParameter(...   
        blk_radiativeHeatTransfer, "area", string(param.facesArea(faceId)),...
        blk_radiativeHeatTransfer, "rad_tr_coeff", string(param.krad(faceId)),...
        blk_thermalMassBlanket, "mass", string(param.facesBlanketMass(faceId)),...
        blk_thermalMassBlanket, "sp_heat", string(param.facesBlanketCp(faceId)));
end

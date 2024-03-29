%% Add all folders to path
% addpath("PlotFunctions/")
% addpath("AuxFunctions/")
% addpath("OrbitPropagators/")
% addpath("Textures/")
% addpath("ThermalModelsSimscape/")


%% ########## Orbit Tab ########
app.param.orb.centralBody = "Earth";
app.param.orb.solarFlux = 1361; % W/m2
app.param.orb.sampleTime = 60; % (s)
app.param.orb.sat = [];
app.param.orb.startDate = datetime(date);
app.param.orb.startHour = 0;
app.param.orb.startMinute = 0;
app.param.orb.durationHours = 1;
app.param.orb.semiMajorAxis = 6778.137; % km
app.param.orb.eccentricity = 0;
app.param.orb.inclination = 52.580; % deg
app.param.orb.RAAN = 190.158; % deg
app.param.orb.argumentOfPeriapsis = 0; % deg
app.param.orb.trueAnomaly = 32.690; % deg
app.param.orb.prop.flag = 0;    
app.param.Tir_Day = 254;% default: 254; % Black-body temp during day time (K)
app.param.Tir_Night = 254;% default: 254; % Black-body temp during night time (K)
app.param.Albedo = 0.3;
app.param.earthIRmodelList = ["Day/Night Temp"];
app.param.moonIRmodelList = ["Day/Night Temp"; "Gradient"];
app.param.IRmodelItemList = app.param.earthIRmodelList;
app.param.IRmodel = "Day/Night Temp";
app.param.moonOrbPropagatorList = ["Kepler"; "Numerical (high precision)"];
app.param.earthOrbPropagatorList = ["two-body-keplerian"; "sgp4"];
app.param.propagatorItemList = app.param.earthOrbPropagatorList;
app.param.orb.orbitPropagator = "two-body-keplerian";
app.param.plot.orbitPlot.showSunVectorFlag = 0;

% Update fields value with inital parameters
app.CentralBodyDropDown.Items = ["Earth"; "Moon"];
app.IRmodelDropDown.Items = app.param.earthIRmodelList;
app.CentralBodyDropDown.Value = app.param.orb.centralBody;
app.IRmodelDropDown.Value = app.param.IRmodel;
app.EpochDatePicker.Value = app.param.orb.startDate;
app.SemimajoraxiskmEditField.Value = app.param.orb.semiMajorAxis;
app.EccentricityEditField.Value = app.param.orb.eccentricity;
app.InclinationEditField.Value = app.param.orb.inclination;
app.RAANEditField.Value = app.param.orb.RAAN;
app.ArgumentofPeriapsisEditField.Value = app.param.orb.argumentOfPeriapsis;
app.TrueAnomalyEditField.Value = app.param.orb.trueAnomaly;
app.SolarfluxWm2EditField.Value = app.param.orb.solarFlux;
app.DayTempKEditField.Value = app.param.Tir_Day;
app.NightTempKEditField.Value = app.param.Tir_Night;
app.AlbedoEditField.Value = app.param.Albedo;
app.PropagatorDropDown.Items = app.param.earthOrbPropagatorList;
app.PropagatorDropDown.Value = app.param.orb.orbitPropagator;
app.UpdateOrbitStatusLabel.Visible = "off";

% Tooltips
packingFactorToolTip = "Percentage of the total panel area that is covered by solar cells";


%% ########## Model Tab ########
app.param.numberOfUnits = 1;

itemList = getPropertiesFromMaterial("options");
app.param.facesMaterial = [itemList(1), itemList(1), itemList(1), itemList(1), itemList(1), itemList(1)];
app.param.satelliteMass = 1;
app.param.solarCellEff = [.28, .28, .28, .28, .28, .28];
app.param.effectiveAreaRatio = [.6, .6, .6, .6, .6, .6];
app.param.conductionResistance1U = 20;   % K/W That is a estimated of the thermal resistance between the internal node and each 1U face
app.param.resistance = ones(1,6)*app.param.conductionResistance1U; % K/W 
app.param.internalNodeCp = 903; %J/kg/K
app.param.internalNode2TotalMassRatio = .7;
app.param.facesMassDistribution = [1/6, 1/6, 1/6, 1/6, 1/6, 1/6]; 
app.param.facesArea = [1/6, 1/6, 1/6, 1/6, 1/6, 1/6]*(6*(0.1*0.1));
app.param.attitude.nadirFace = "X+";
app.param.attitude.ramFace = "Y+";
app.param.attitude.ramFaceItems = ["Y+", "Y-", "Z+", "Z-"];
app.param.deployable.type = ["None"; "None"; "None"; "None"; "None"; "None"];
app.param.deployable.hinge = ["X+ | Y+"; "X- | Y+"; "Y+ | X+"; "Y- | X+"; "Z+ | X+"; "Z- | X+"];
app.param.deployable.size = ["1P"; "1P"; "1P"; "1P"; "1P"; "1P"];
app.param.deployable.solarCellEff = [.28, .28, .28, .28, .28, .28];
app.param.deployable.effectiveAreaRatio = [.6, .6, .6, .6, .6, .6];
app.param.deployable.flipFixedPanel = [0,0,0,0,0,0];

% Mass Table
massVector = app.param.satelliteMass * [app.param.internalNode2TotalMassRatio; (1-app.param.internalNode2TotalMassRatio)* app.param.facesMassDistribution'];
app.UITableMassDistribution.Data = table(["Internal Node";"X+";"X-";"Y+";"Y-";"Z+";"Z-"],massVector,'VariableNames',{'Node','Mass'});


% Insulator (Blanket) Properties
app.param.resistanceBlanket = [0, 0, 0, 0, 0, 0];                   % Thermal resistance between face inner surface and the outer surface of the blanket K/W 
app.param.facesBlanketMass = 1e-5 * ones(1,6);                      % Blanket mass at each face (kg)
app.param.facesBlanketCp = 130 * ones(1,6);                         % Gold specific heat to the blanket

% The emissivity for the heat transfer though MLI is computed considering
% the emissivity of the face surface (Aluminum - Rough) and the MLI
e_Al = 0.07; e_MLI = 0.025;
emissivity_Al2MLI = 1/(1/e_Al + 1/e_MLI - 1); % Ref.: Eq. (3) in https://www.sciencedirect.com/science/article/pii/S092145342030397X
sigma = 5.67e-8; % Stefan-Boltzmann constant (W · m -2 · K -4 )
app.param.kradFace2Blanket = ones(1,6)*sigma*emissivity_Al2MLI; 
app.param.resistanceInterFacesBlanket = ones(6,6) * 1e3;                   % Matrix of inter-faces (blanket) resistance K/W 

% Face properties
app.param.absorptivity = [getPropertiesFromMaterial(app.param.facesMaterial(1)).absorp,...
                        getPropertiesFromMaterial(app.param.facesMaterial(2)).absorp,...
                        getPropertiesFromMaterial(app.param.facesMaterial(3)).absorp,...
                        getPropertiesFromMaterial(app.param.facesMaterial(4)).absorp,...
                        getPropertiesFromMaterial(app.param.facesMaterial(5)).absorp,...
                        getPropertiesFromMaterial(app.param.facesMaterial(6)).absorp]; 
app.param.emissivity = [getPropertiesFromMaterial(app.param.facesMaterial(1)).emiss,...
                        getPropertiesFromMaterial(app.param.facesMaterial(2)).emiss,...
                        getPropertiesFromMaterial(app.param.facesMaterial(3)).emiss,...
                        getPropertiesFromMaterial(app.param.facesMaterial(4)).emiss,...
                        getPropertiesFromMaterial(app.param.facesMaterial(5)).emiss,...
                        getPropertiesFromMaterial(app.param.facesMaterial(6)).emiss];
app.param.facesCp = [getPropertiesFromMaterial(app.param.facesMaterial(1)).Cp,...
                        getPropertiesFromMaterial(app.param.facesMaterial(2)).Cp,...
                        getPropertiesFromMaterial(app.param.facesMaterial(3)).Cp,...
                        getPropertiesFromMaterial(app.param.facesMaterial(4)).Cp,...
                        getPropertiesFromMaterial(app.param.facesMaterial(5)).Cp,...
                        getPropertiesFromMaterial(app.param.facesMaterial(6)).Cp]; 
app.param.resistanceInterFaces = ones(6,6) * 1e3;                   % Matrix of inter-faces (chassis) resistance K/W 

% Update X+ face tab
app.SurfaceTypeDropDownXplus.Items = itemList;
p = getPropertiesFromMaterial(itemList(1));
app.EmissivitySpinnerXplus.Value = p.emiss;
app.EmissivitySpinnerXplus.Enable = "off";
app.AbsorptivitySpinnerXplus.Value = p.absorp;
app.AbsorptivitySpinnerXplus.Enable = "off";
app.RatioEditFieldXplus.Value = p.absorp/p.emiss;
app.CpSpinnerXplus.Value = p.Cp;
app.CpSpinnerXplus.Enable = "off";
app.HingeDropDownXplus.Enable = "off";
app.HingeDropDownXplus.Items = ["X+ | Y+"; "X+ | Y-"; "X+ | Z+"; "X+ | Z-"];
app.SizeDropDownXplus.Items = ["1P", "2P", "3P", "4P","5P","6P"];
app.SizeDropDownXplus.Enable = "off";
app.CellEfficiencySpinnerXplusDeployable.Enable = "off";
app.EffecAreaRatioSpinnerXplusDeployable.Enable = "off";
app.FlipPanelCheckBoxXplus.Enable = "off";
app.EffecAreaRatioSpinnerXplus.Tooltip = packingFactorToolTip;
app.EffecAreaRatioSpinnerXplusDeployable.Tooltip = packingFactorToolTip;


% Update X- face tab
app.SurfaceTypeDropDownXminus.Items = itemList;
p = getPropertiesFromMaterial(itemList(1));
app.EmissivitySpinnerXminus.Value = p.emiss;
app.EmissivitySpinnerXminus.Enable = "off";
app.AbsorptivitySpinnerXminus.Value = p.absorp;
app.AbsorptivitySpinnerXminus.Enable = "off";
app.RatioEditFieldXminus.Value = p.absorp/p.emiss;
app.CpSpinnerXminus.Value = p.Cp;
app.CpSpinnerXminus.Enable = "off";
app.HingeDropDownXminus.Enable = "off";
app.HingeDropDownXminus.Items = ["X- | Y+"; "X- | Y-"; "X- | Z+"; "X- | Z-"];
app.SizeDropDownXminus.Items = ["1P", "2P", "3P", "4P","5P","6P"];
app.SizeDropDownXminus.Enable = "off";
app.CellEfficiencySpinnerXminusDeployable.Enable = "off";
app.EffecAreaRatioSpinnerXminusDeployable.Enable = "off";
app.FlipPanelCheckBoxXminus.Enable = "off";
app.EffecAreaRatioSpinnerXminus.Tooltip = packingFactorToolTip;
app.EffecAreaRatioSpinnerXminusDeployable.Tooltip = packingFactorToolTip;

% Update Y+ face tab
app.SurfaceTypeDropDownYplus.Items = itemList;
p = getPropertiesFromMaterial(itemList(1));
app.EmissivitySpinnerYplus.Value = p.emiss;
app.EmissivitySpinnerYplus.Enable = "off";
app.AbsorptivitySpinnerYplus.Value = p.absorp;
app.AbsorptivitySpinnerYplus.Enable = "off";
app.RatioEditFieldYplus.Value = p.absorp/p.emiss;
app.CpSpinnerYplus.Value = p.Cp;
app.CpSpinnerYplus.Enable = "off";
app.HingeDropDownYplus.Enable = "off";
app.HingeDropDownYplus.Items = ["Y+ | X+"; "Y+ | X-"; "Y+ | Z+"; "Y+ | Z-"];
app.SizeDropDownYplus.Items = ["1P", "2P", "3P", "4P","5P","6P"];
app.SizeDropDownYplus.Enable = "off";
app.CellEfficiencySpinnerYplusDeployable.Enable = "off";
app.EffecAreaRatioSpinnerYplusDeployable.Enable = "off";
app.FlipPanelCheckBoxYplus.Enable = "off";
app.EffecAreaRatioSpinnerYplus.Tooltip = packingFactorToolTip;
app.EffecAreaRatioSpinnerYplusDeployable.Tooltip = packingFactorToolTip;

% Update Y- face tab
app.SurfaceTypeDropDownYminus.Items = itemList;
p = getPropertiesFromMaterial(itemList(1));
app.EmissivitySpinnerYminus.Value = p.emiss;
app.EmissivitySpinnerYminus.Enable = "off";
app.AbsorptivitySpinnerYminus.Value = p.absorp;
app.AbsorptivitySpinnerYminus.Enable = "off";
app.RatioEditFieldYminus.Value = p.absorp/p.emiss;
app.CpSpinnerYminus.Value = p.Cp;
app.CpSpinnerYminus.Enable = "off";
app.HingeDropDownYminus.Enable = "off";
app.HingeDropDownYminus.Items = ["Y- | X+"; "Y- | X-"; "Y- | Z+"; "Y- | Z-"];
app.SizeDropDownYminus.Items = ["1P", "2P", "3P", "4P","5P","6P"];
app.SizeDropDownYminus.Enable = "off";
app.CellEfficiencySpinnerYminusDeployable.Enable = "off";
app.EffecAreaRatioSpinnerYminusDeployable.Enable = "off";
app.FlipPanelCheckBoxYminus.Enable = "off";
app.EffecAreaRatioSpinnerYminus.Tooltip = packingFactorToolTip;
app.EffecAreaRatioSpinnerYminusDeployable.Tooltip  = packingFactorToolTip;

% Update Z+ face tab
app.SurfaceTypeDropDownZplus.Items = itemList;
p = getPropertiesFromMaterial(itemList(1));
app.EmissivitySpinnerZplus.Value = p.emiss;
app.EmissivitySpinnerZplus.Enable = "off";
app.AbsorptivitySpinnerZplus.Value = p.absorp;
app.AbsorptivitySpinnerZplus.Enable = "off";
app.RatioEditFieldZplus.Value = p.absorp/p.emiss;
app.CpSpinnerZplus.Value = p.Cp;
app.CpSpinnerZplus.Enable = "off";
app.HingeDropDownZplus.Enable = "off";
app.HingeDropDownZplus.Items = ["Z+ | X+"; "Z+ | X-"; "Z+ | Y+"; "Z+ | Y-"];
app.SizeDropDownZplus.Items = ["1P", "2P", "3P", "4P","5P","6P"];
app.SizeDropDownZplus.Enable = "off";
app.CellEfficiencySpinnerZplusDeployable.Enable = "off";
app.EffecAreaRatioSpinnerZplusDeployable.Enable = "off";
app.FlipPanelCheckBoxZplus.Enable = "off";
app.EffecAreaRatioSpinnerZplus.Tooltip = packingFactorToolTip;
app.EffecAreaRatioSpinnerZplusDeployable.Tooltip = packingFactorToolTip;

% Update Z- face tab
app.SurfaceTypeDropDownZminus.Items = itemList;
p = getPropertiesFromMaterial(itemList(1));
app.EmissivitySpinnerZminus.Value = p.emiss;
app.EmissivitySpinnerZminus.Enable = "off";
app.AbsorptivitySpinnerZminus.Value = p.absorp;
app.AbsorptivitySpinnerZminus.Enable = "off";
app.RatioEditFieldZminus.Value = p.absorp/p.emiss;
app.CpSpinnerZminus.Value = p.Cp;
app.CpSpinnerZminus.Enable = "off";
app.HingeDropDownZminus.Enable = "off";
app.HingeDropDownZminus.Items = ["Z- | X+"; "Z- | X-"; "Z- | Y+"; "Z- | Y-"];
app.SizeDropDownZminus.Items = ["1P", "2P", "3P", "4P","5P","6P"];
app.SizeDropDownZminus.Enable = "off";
app.CellEfficiencySpinnerZminusDeployable.Enable = "off";
app.EffecAreaRatioSpinnerZminusDeployable.Enable = "off";
app.FlipPanelCheckBoxZminus.Enable = "off";
app.EffecAreaRatioSpinnerZminus.Tooltip = packingFactorToolTip;
app.EffecAreaRatioSpinnerZminusDeployable.Tooltip = packingFactorToolTip;


% Attitude
app.NadirFaceDropDown.Items = ["X+"; "X-"; "Y+"; "Y-"; "Z+"; "Z-"];
app.NadirFaceDropDown.Value = app.param.attitude.nadirFace;
app.RamFaceDropDown.Items = ["Y+"; "Y-"; "Z+"; "Z-"];
app.RamFaceDropDown.Value = app.param.attitude.ramFace;
app.RamFaceDropDown.Tooltip = "Ram side is the side that points in the direction of the satellite's motion";

% Internal node properties
app.IntermalnodetototalmassratioSpinner.Tooltip = "Ratio of the internal nodel mass to the total satellite mass. The remaining mass is distributed to face nodes proportionally to the each face area.";

% Thermal Resistance
app.SpinnerRXplus.Tooltip = "Thermal resistance between the internal nodel and the X+ face";
app.SpinnerRXminus.Tooltip = "Thermal resistance between the internal nodel and the X- face";
app.SpinnerRYplus.Tooltip = "Thermal resistance between the internal nodel and the Y+ face";
app.SpinnerRYminus.Tooltip = "Thermal resistance between the internal nodel and the Y- face";
app.SpinnerRZplus.Tooltip = "Thermal resistance between the internal nodel and the Z+ face";
app.SpinnerRZminus.Tooltip = "Thermal resistance between the internal nodel and the Z- face";
app.SpinnerRXplus.Value = app.param.resistance(1);
app.SpinnerRXminus.Value = app.param.resistance(2);
app.SpinnerRYplus.Value = app.param.resistance(3);
app.SpinnerRYminus.Value = app.param.resistance(4);
app.SpinnerRZplus.Value = app.param.resistance(5);
app.SpinnerRZminus.Value = app.param.resistance(6);
app.ThermalResistanceImage.ImageSource = "Figures/ThermalModel7NodesSchematic_v2.png";


%% ########## Power Tab ########

app.param.pwrDissipationProfileList = ["Constant"; "Day/Night"];
app.param.pwrDissipationProfile = app.param.pwrDissipationProfileList(1);
app.param.pwrDissipationDaySimpleModel = 1;      % W
app.param.pwrDissipationNightSimpleModel = 1;    % W
app.param.pwr.flag = 0;
app.param.computePwrWithShadowFlag = 0;
app.param.plot.animation.status = "stop";

% Update fields value with inital parameters
app.ProfileDropDown.Items = app.param.pwrDissipationProfileList;
app.ProfileDropDown.Value = app.param.pwrDissipationProfile;
app.PwrDayWSpinner.Value = app.param.pwrDissipationDaySimpleModel;
app.PwrNightWSpinner.Value = app.param.pwrDissipationNightSimpleModel;
app.PwrNightWSpinner.Enable = "off";
app.ComputePowerStatusLabel.Visible = "off";
app.ConsidershadowsCheckBox.Value = app.param.computePwrWithShadowFlag;

%% ########## Thermal Tab ########
app.param.thermal.flag = 0; 
app.param.plot.solarRadiationInput.showFacesFlag = 0;
app.param.plot.IRRadiationInput.showFacesFlag = 0;
app.param.plot.albedoRadiationInput.showFacesFlag = 0;
app.param.plot.IRRadiationInput.showAbsorbedIRFlag = 1;
app.param.plot.IRRadiationInput.showRadiatedIRFlag = 0;

% Update fields value with inital parameters
app.ComputeThermalStatusLabel.Visible = "off";



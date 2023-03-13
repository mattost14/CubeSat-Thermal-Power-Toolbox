function p = getPropertiesFromMaterial(material)
    switch material
        case 'options'
            p = ["Gold coating"; "Black coating"; "White coating"; "Solar Panel"; "Custom"];
        case 'Solar Panel'
            p.Cp = 1164; % Specific heat (J/(kg*K)) FR4: 1164 J/(kg*K)
            p.rho = 1809; % Density (kg/m3)
            p.absorp = 0.92*.6+.76*.4; %  Solar Panel Absorptivity = 60%(Solar Cell Absorptivity) + 40%(Al-6061-T6-Black Anodized/FR4 Absorptivity)
            p.emiss = 0.85*.6+0.88*.4; % Emissivity Solar
        case 'Black coating'
            p.Cp = 960; % Specific heat (J/(kg*K))
            p.rho = 2765; % Density (kg/m3)
            p.absorp = 0.96; % Absorptivity Solar
            p.emiss = 0.91; % Emissivity Solar
        case 'White coating'
            p.Cp = 960; % Specific heat (J/(kg*K))
            p.rho = 2765; % Density (kg/m3)
            p.absorp = 0.2; % Absorptivity Solar
            p.emiss = 0.92; % Emissivity Solar
        case 'Gold coating'
            p.Cp = 960; % Specific heat (J/(kg*K))
            p.rho = 2765; % Density (kg/m3)
            p.absorp = 0.19; % Absorptivity Solar
            p.emiss = 0.025; % Emissivity Solar 
    end
end
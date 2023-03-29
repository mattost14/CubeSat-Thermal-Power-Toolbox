function p = getPropertiesFromMaterial(material)
    switch material
        case 'options'
            p = ['Aluminum (Polished)'; "Black coating"; "White coating"; "Solar Panel"; "MLI blanket";"Custom"];
        case 'blankets'
            p = ["MLI blanket"];
        case 'Solar Panel'
            p.Cp = 1164; % Specific heat (J/(kg*K)) FR4: 1164 J/(kg*K)
            p.areaDensity = 1809; % Density (kg/m3)
            p.absorp = 0.92*.6+.76*.4; %  Solar Panel Absorptivity = 60%(Solar Cell Absorptivity) + 40%(Al-6061-T6-Black Anodized/FR4 Absorptivity)
            p.emiss = 0.85*.6+0.88*.4; % Emissivity Solar
        case 'Black coating'
            p.Cp = 960; % Specific heat (J/(kg*K))
            p.absorp = 0.96; % Absorptivity Solar
            p.emiss = 0.91; % Emissivity Solar
        case 'White coating'
            p.Cp = 960; % Specific heat (J/(kg*K))
            p.absorp = 0.2; % Absorptivity Solar
            p.emiss = 0.92; % Emissivity Solar
        case 'MLI blanket'
            p.Cp = 960; % Specific heat (J/(kg*K))
            p.areaDensity = 1.5; % Area Density (kg/m2)
            p.absorp = 0.19; % Absorptivity Solar
            p.emiss = 0.025; % Emissivity Solar
            p.conductionResistanceTimesArea = 400; % Conduction Resistance to cross 1m2 of MLI (K*m^2/W) - That is the estimated thermal resistance to cross 10 layers of MLI in 1U face area (Ref: https://www2.jpl.nasa.gov/adv_tech/coolers/Cool_ppr/CEC2015%20Quantifying%20MLI%20Thermal%20Conduction.pdf)
        case 'Aluminum (Polished)'
            p.Cp = 960; % Specific heat (J/(kg*K))
            p.absorp = 0.15; % Absorptivity Solar
            p.emiss = 0.05; % Emissivity Solar
    end
end
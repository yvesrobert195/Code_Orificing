function G=geometry(Core)
if strcmpi(Core,'A_one_row')
    Core='A';
end
if strcmpi(Core,'A')
    %% Geometry [m],[Pa]
    % Rods geometry
    G.Rod.pitch=1.3608e-2; % Space between between rods
    G.Rod.fuel_diameter=0.9516e-2;
    G.Rod.rod_diameter=1.2208e-2;
    G.Rod.clad_thickness=0.6104e-3;
    G.Rod.height=3.18;
    
    % Assemblies geometry
    G.Assembly.pins_per_assembly=169;
    G.Assembly.duct_thickness=3.5e-3;
    G.Assembly.gap=2e-3; % Inter-assembly gap
    G.Assembly.pitch=22.05e-2; % Space between between assemblies
    G.Assembly.flow_area=123.3184E-4; % Flow area per assembly, m^2
    G.Assembly.loss_coefficient=1.29+2.01+0.363+0.41+0.098+0.37+0.79+1.32; % Pressure loss coefficient for each assembly
    G.Assembly.wetted_perimeter = 7.2062; %wetted perimeter per assembly, m
    G.Assembly.hydraulic_diameter=4*G.Assembly.flow_area/G.Assembly.wetted_perimeter;
    
    % Core geometry
    G.Core.number_assemblies=196;
    G.Core.fuel_pins=33124;
    G.Core.active_volume=18.58;
    G.Core.equivalent_diameter=280.79e-2;
    G.Core.active_height=300e-2;
elseif strcmpi(Core,'bnb')
    %% Geometry [m],[Pa]
    % Rods geometry
    G.Rod.pitch=1.3609e-2; % Space between between rods
    G.Rod.fuel_diameter=0.9523e-2;
    G.Rod.rod_diameter=1.2217e-2;
    G.Rod.clad_thickness=0.6109e-3;
    G.Rod.height=3.18;
    
    % Assemblies geometry
    G.Assembly.pins_per_assembly=169;
    G.Assembly.duct_thickness=3.5e-3;
    G.Assembly.gap=2e-3; % Inter-assembly gap
    G.Assembly.pitch=22.05e-2; % Space between between assemblies
    G.Assembly.flow_area=123.3184E-4; % Flow area per assembly, m^2
    G.Assembly.loss_coefficient=1.29+2.01+0.363+0.41+0.098+0.37+0.79+1.32; % Pressure loss coefficient for each assembly
    G.Assembly.wetted_perimeter = 7.2062; %wetted perimeter per assembly, m
    G.Assembly.hydraulic_diameter=4*G.Assembly.flow_area/G.Assembly.wetted_perimeter;
    
    % Core geometry
    G.Core.number_assemblies=492;
    G.Core.fuel_pins=83148;
    G.Core.active_volume=46.63;
    G.Core.equivalent_diameter=444.87e-2;
    G.Core.active_height=300e-2;
else
    error('Please specify a valid core name')
end
% %% TADPOLE Torch Sizing
% % Authors: Jonah Nobel
% % First Created: 11/11/2023
% % Last Updated: 11/11/2023
% % Calculations done in SI units

%%% Definitions

%%Inputs
chamber_pressure_psi = 150;                                                 %Psia
exit_pressure_psi = 14.7;                           
ox_feed_pressure_psi = 550;
fuel_feed_pressure_psi = 550;
percent_chamber_m_dot = 1.5;                                                % % of mass flow
fuel_inlet_diameter_in = 1/4;
c_d = 0.7;                                                                  %Discharge coeficcient from literature (Sutton, p279)
OF = 1.0;

%%Selected orifices
%Possible orifice sizes: https://www.mcmaster.com/products/orifices/threaded-flow-control-orifices/body-material~stainless-steel-1/body-material~316-stainless-steel/thread-type~unf-1/
chosen_ox_diameter_in = 0.018;                                              %in
chosen_fuel_diameter_in = 0.020;                                            %in

fuel = "C3H8O,2propanol";
oxidizer = "O2";


%Pressures and OF
P_c = chamber_pressure_psi * 6894.757;                                      %Psi to Pascals
P_e = 14.7 * 6894.757;                              
ox_feed_pressure = ox_feed_pressure_psi * 6894.757; 
fuel_feed_pressure = fuel_feed_pressure_psi * 6894.757;
fuel_inlet_area = ((fuel_inlet_diameter_in/2) * 0.5 * 0.0254)^2 * pi;       %m^2

%Constants
CEA_input_name = 'torchfile';
R = 8.314;

%Propellant information
fuel_weight = 0;
fuel_density = 786;                                                         %kg/m^3
fuel_temp = 293.15;                                                         %Kelvin
fuel_gamma = 154.75 / (154.75-R);                                           %https://webbook.nist.gov/cgi/cbook.cgi?ID=C67630&Mask=7#Thermo-Condensed

oxidizer_temp = 293.15;                                                     %Kelvin
ox_density = 50.1;                                                          %@ 550 psi, change for actual feed pressure
ox_gamma = 1.40;

%m_dot_definitions
m_dot_chamber = 2.72 * 0.4535924;                                           %kg
m_dot_torch = (percent_chamber_m_dot/100) * m_dot_chamber;                                    %m_dot_chamber * 0.015  * 0.4535924 %Kilograms
m_dot_fuel = m_dot_torch/(1+OF);                                            %kg
m_dot_ox = m_dot_fuel*OF;                                                   %kg

%Chosen orifice diameters
chosen_ox_diameter = chosen_ox_diameter_in*0.0254;                          %in to m
chosen_fuel_diameter = chosen_fuel_diameter_in*0.0254;                      %in to m

%%% Calculations

[c_star_t, ~, ~, M_t, gamma_t, P_g_t, T_g_t, ~, mu_g_t, Pr_g_t, ~, ~, ~, cp_g_t] = RunCEA(P_c, P_e, fuel, fuel_weight, fuel_temp, oxidizer, oxidizer_temp, OF, 0, 0, 0, 0, 0, CEA_input_name);

A_t = c_star_t * m_dot_torch / P_c;                                         %In SI units
A_t_inches = A_t * 1550;                                                    %inches^2
throat_diameter = sqrt(A_t_inches / pi) * 2;                                %inches

ox_choked = (2/(ox_gamma + 1))^(ox_gamma/(ox_gamma-1)) > (P_c/ox_feed_pressure);

ox_area = (m_dot_ox * sqrt(oxidizer_temp)) / (ox_feed_pressure * sqrt(ox_gamma/R) * ((ox_gamma+1)/2)^(-1*(ox_gamma+1)/(2*(ox_gamma-1)))); %m^2, NASA Choked flow equation, M = 1
ox_area_inches = ox_area * 1550;                                            %inches^2
ox_diameter_inches = sqrt(ox_area_inches / pi) * 2;                         %inches

fuel_area = m_dot_fuel/(c_d*sqrt(2*fuel_density*(fuel_feed_pressure-P_c))); %m^2, orifice sizing equation
fuel_area_inches = fuel_area * 1550;                                        %inches^2
fuel_diameter_inches = sqrt(fuel_area_inches / pi) * 2;                     %inches  


%Calculations based on chosen orifices
new_ox_area = (chosen_ox_diameter/2)^2 * pi;
new_fuel_area = (chosen_fuel_diameter/2)^2 * pi;
new_ox_m_dot = (((new_ox_area) * ox_feed_pressure)/sqrt(oxidizer_temp)) * sqrt(ox_gamma/R) * ((ox_gamma+1)/2)^(-1*(ox_gamma+1)/(2*(ox_gamma-1))); %Assume choked flow
new_fuel_m_dot = c_d  * new_fuel_area * sqrt(2*fuel_density*(fuel_feed_pressure-P_c));


new_of = new_ox_m_dot / new_fuel_m_dot;
[new_c_star, ~, ~, new_M_t, new_gamma_t, new_P_g_t, new_T_g_t, ~, new_mu_g_t, new_Pr_g_t, ~, ~, ~, new_cp_g_t] = RunCEA(P_c, P_e, fuel, fuel_weight, fuel_temp, oxidizer, oxidizer_temp, new_of, 0, 0, 0, 0, 0, CEA_input_name);

new_A_t = new_c_star * (new_ox_m_dot+ new_fuel_m_dot) / P_c;                %In SI units
new_A_t_inches = new_A_t * 1550;                                            %inches^2
new_throat_diameter = sqrt(new_A_t_inches / pi) * 2;                        %inches


fprintf("\nCalculations from chamber pressure:\n")
fprintf("Chamber throat diameter is %0.4f inches.\n", throat_diameter)
if ox_choked == 1
    fprintf("Ox is choked.\n")
elseif ox_choked == 0
   fprintf("Ox is unchoked.\n")
end
fprintf("Oxidizer injection orifice area is %f in^2, and ideal oxidizer orifice diameter is %0.4f inches.\n", ox_area_inches, ox_diameter_inches)
fprintf("Fuel injection orifice area is %f in^2, and ideal fuel orifice diameter is %0.4f inches.\n", fuel_area_inches, fuel_diameter_inches)

fprintf("\nSelected orifice diameters are %.4f in for ox and %.4f in for fuel.\n", chosen_ox_diameter_in, chosen_fuel_diameter_in)
fprintf("Calculations from chosen orifice diameters and feed pressure:\n")
fprintf("Ox massflow is %.4f kg/s, and fuel mass flow is %.4f kg/s.\n", new_ox_m_dot, new_fuel_m_dot)
fprintf("Torch mass flow is %.4f kg/s, %.2f%% of total engine mass flow.\n", new_ox_m_dot + new_fuel_m_dot, 100*(new_ox_m_dot + new_fuel_m_dot)/m_dot_chamber)
fprintf("OF ratio is now %.4f at c-star %.0f.\n", new_of, new_c_star)
fprintf("New throat diameter is %.4f inches at chamber pressure %.1f psia.\n", new_throat_diameter, chamber_pressure_psi) 
%% Thermal  M odel for Nano Satellite
%% Physical constants
clc; clear all;
qe = 240; % radian energy from object during eclipse
qs =1400; % solar radiation constant
qa =0.33; % albedo radiation constant
qp = 245; % planetary radiation constant
sigma = 5.67e-8; % Stefan-Boltzmann constant [W/m^2.K^4]
RE = 6371e3; % Earth radius [m]
Pdis= 0.2; % power dissipation
%% Satellite specifications
mass = 1; % [kg]
alpha = 0.3; % Absorptance
epsilon = 0.08; % Emmitance
epsilon_earth = 0.6; % Emmitance of earth
Ti = 15 + 273.15; % initial temperature [k]
cp = 0.902; % Specific heat capacity [J/kgk]
h = 408e3; % Altitude [m]
a = sqrt(10e-4); % Frontal area is 10 cm^2, finding side length of cube
A = 6*a^2; % total area [m^2]
As = a^2; % Frontal area [m^2]
Ap = a^2;
Aa = a^2;
%% Thermal time constant
tau_t =(mass * cp) / (epsilon * sigma * A);
%% The eclipse and orbital duration
theta = asind(RE / (RE + h));
P = 2 * pi * sqrt((RE + h)^3 / (6.67430e-11 * 5.972e24)); % Orbital period using Kepler's third law
te = P * (theta / 360); % Eclipse time
ne = P-te;
F12 = cos(theta)/(1+RE+h)^2; % view factor

% Thermal Power Input during sun exposure
phi = 0;
theta = 0;
Qin_sun = qs*alpha*As + qa*qs*alpha*F12*Ap*cos(phi)*cos(theta) + qp*epsilon_earth*F12*Ap+Pdis ;

% Equilibrium Temperature during sun exposure
Teq_sun = (Qin_sun / (epsilon * sigma * A))^(1/4);

% During eclipse
Qin_eclipse = qe * A * epsilon + Pdis ; % No solar radiation during eclipse
Teq_eclipse = (Qin_eclipse / (epsilon * sigma * A))^(1/4);
dt = 1; % Time step, e.g., 1 second
t = te; % Starting from the end of the eclipse period
T = Teq_eclipse; % Starting temperature at the end of the eclipse
T_array = [];
time_array = [];
counter = 1;
%% Three orbital rotations
while t < te + (3*P)
    if t >= te && t < te + P % First sunlight period after eclipse
        dT = ((Teq_sun^4 - T^4) / tau_t) * dt;
    elseif t >= te + P && t < te + 2*P % Second orbital period (including eclipse)
        if t < te + P + te % If within the second eclipse duration
            dT = ((Teq_eclipse^4 - T^4) / tau_t) * dt;
        else
            dT =( (Teq_sun^4 - T^4) / tau_t )* dt;
        end
    else % Third orbital period (including eclipse)
        if t < te + 2*P + te % If within the third eclipse duration
            dT = (Teq_eclipse^4 - T^4) / tau_t * dt;
        else
            dT = (Teq_sun^4 - T^4) / tau_t * dt;
        end
    end

    T = T + dT;
    t = t + dt;

    T_array(counter) = T;
    time_array(counter) = t - te; % Adjusting time to start from 0
    counter = counter + 1;
end
%% Plot
plot(time_array/60, T_array - 273.15, 'black-'); % Convert seconds to minutes for plotting
xlabel('Time [minutes]');
ylabel('Temperature [°C]');
title('Satellite Temperature Cycle');
grid on;

% Highlight eclipse periods
hold on;
% Start shading from the second orbit
for i = 1:3
    x_fill = [i*P/60, i*P/60, i*P/60 + te/60, i*P/60 + te/60];
    y_fill = [min(T_array-273.15), 140, 140, min(T_array-273.15)];
    fill(x_fill, y_fill, [0.9,0.9,0.9], 'FaceAlpha', 1, 'EdgeColor','none');
end
legend('Temperature', 'Eclipse');
% If the eclipse regions are covering the temperature curve, send them to the back.
uistack(findall(gca, 'Type', 'patch'), 'bottom');
hold off;

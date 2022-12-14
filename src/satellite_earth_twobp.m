%% Two-body problem ODE
function state_dot = satellite_earth_twobp(~, state)
    ME = 5.97217e+24;   % Mass of earth [kg]
    G =  6.67430e-11;   % Gravitational constant [m³/(kg.s²)]
    GM = G*ME;          % Earth's: gravitational constant * mass [m³/s²]
    x = state(1);
    y = state(2);
    z = state(3);
    x_dot = state(4);
    y_dot = state(5);
    z_dot = state(6);
    x_ddot = -GM * x / (x^2 + y^2 + z^2)^(3/2);
    y_ddot = -GM * y / (x^2 + y^2 + z^2)^(3/2);
    z_ddot = -GM * z / (x^2 + y^2 + z^2)^(3/2);
    state_dot = [x_dot; y_dot; z_dot; x_ddot; y_ddot; z_ddot];
end

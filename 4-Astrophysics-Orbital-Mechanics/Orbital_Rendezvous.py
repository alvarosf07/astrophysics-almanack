global mu;
global r_earth;

mu = 3.986e5;% [km^3/s2]
r_earth = 6378.14; % [km]


% Exercise 1
fprintf ("\n*********************************************\n");
fprintf ("\n Third Exercise - Orbit Phasing & Rendezvous\n");
fprintf ("\n*********************************************\n");

%%
% A) Calculation of the ISS Initial State Vector
fprintf ("\nA&B) Calculation of the ISS Initial State Vector\n");

% Initial Data
delta_theta = 100*pi/180;
chaser_ISS = 180*pi/180;
theta_ISS = chaser_ISS+delta_theta;
r_ISS = r_earth+404;

% Calculation of the ISS Initial State Vector
rx_ISS = -r_ISS*cos(theta_ISS);
ry_ISS = -r_ISS*sin(theta_ISS);
rz_ISS = 0;
r0_ISS = [rx_ISS,ry_ISS,rz_ISS];

v_ISS = sqrt (mu/r_ISS);
vx_ISS = -v_ISS*cos(delta_theta-pi/2);
vy_ISS = -v_ISS*sin(delta_theta-pi/2);
vz_ISS = 0;
v0_ISS = [vx_ISS,vy_ISS,vz_ISS];

X0_ISS = [rx_ISS;ry_ISS;rz_ISS;vx_ISS;vy_ISS;vz_ISS];

T_ISS = 2*pi*sqrt(r_ISS^3/mu);

fprintf ("     - rx             (km)         = %g\n", r0_ISS(1));
fprintf ("     - ry             (km)         = %g\n", r0_ISS(2));
fprintf ("     - rz             (km)         = %g\n", r0_ISS(3));
fprintf ("     - vx             (km/s)       = %g\n", v0_ISS(1));
fprintf ("     - vy             (km/s)       = %g\n", v0_ISS(2));
fprintf ("     - vz             (km/s)       = %g\n", v0_ISS(3));
fprintf ("     - T              (s)          = %g\n", T_ISS);

% Plot the Orbit
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
[t,f] = ode45(@twobody,[0:10:5.429877739056294e+03*13],X0_ISS,options);
t_ISS = t';
f_ISS = f';
r_vector_ISS = [f_ISS(1,:);f_ISS(2,:);f_ISS(3,:)];
r_mod_vector_ISS = sqrt( f_ISS(1,:).^2 + f_ISS(2,:).^2 + f_ISS(3,:).^2 );
v_vector_ISS = [f_ISS(4,:);f_ISS(5,:);f_ISS(6,:)];
v_mod_vector_ISS = sqrt( f_ISS(4,:).^2 + f_ISS(5,:).^2 + f_ISS(6,:).^2 );

figure (1);
hold on;
plot3 (f_ISS(1,:),f_ISS(2,:),f_ISS(3,:),'-b','LineWidth',4);
[earth_x,earth_y,earth_z] = sphere(50);
earth_x = earth_x*r_earth + 0;
earth_y = earth_y*r_earth + 0;
earth_z = earth_z*r_earth + 0;
surf (earth_x,-earth_y,-earth_z);
h = findobj('Type', 'surface');
%image = 'http://upload.wikimedia.org/wikipedia/commons/thumb/c/cd/Land_ocean_ice_2048.jpg/1024px-Land_ocean_ice_2048.jpg';
image = imread('Earth_Texture.jpg');
imshow(image);
set(h,'CData',image, 'FaceColor','texturemap','edgecolor','none')
grid off;
axis on;
axis equal;
xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');
view([1 1 0.75]); % adjust the viewing angle
%zoom(2);
hold off;

%%
% B) Calculation of a & e of Chaser Orbit
fprintf ("\n\nB) Calculation of a & e of Chaser Orbit\n");
N_rev = 12;
a_chaser = ( (1 - (delta_theta/(2*pi*N_rev)) )^(2/3) ) * r_ISS;
e_chaser = (r_ISS/a_chaser) - 1.0 ;
rp_chaser = a_chaser*(1-e_chaser);
ra_chaser = a_chaser*(1+e_chaser);

% Calculation of the chaser Initial State Vector
rx_chaser = ra_chaser;
ry_chaser = 0;
rz_chaser = 0;
r0_chaser = [rx_chaser,ry_chaser,rz_chaser];

v_chaser = sqrt ( ((2*mu)/ra_chaser)-(mu/a_chaser) );
vx_chaser = 0;
vy_chaser = v_chaser;
vz_chaser = 0;
v0_chaser = [vx_chaser,vy_chaser,vz_chaser];

X0_chaser = [rx_chaser;ry_chaser;rz_chaser;vx_chaser;vy_chaser;vz_chaser]

T_chaser = 2*pi*sqrt(a_chaser^3/mu)

% Plot the Orbit
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
[t,f] = ode45(@twobody,[0:10:T_chaser*13],X0_chaser,options);
t_chaser = t';
f_chaser = f';
r_vector_chaser = [f_chaser(1,:);f_chaser(2,:);f_chaser(3,:)];
r_mod_vector_chaser = sqrt( f_chaser(1,:).^2 + f_chaser(2,:).^2 + f_chaser(3,:).^2 );
v_vector_chaser = [f_chaser(4,:);f_chaser(5,:);f_chaser(6,:)];
v_mod_vector_chaser = sqrt( f_chaser(4,:).^2 + f_chaser(5,:).^2 + f_chaser(6,:).^2 );

figure (2);
hold on;
plot3 (f_ISS(1,:),f_ISS(2,:),f_ISS(3,:),'-b','LineWidth',2);
plot3 (f_chaser(1,:),f_chaser(2,:),f_chaser(3,:),'-g','LineWidth',2);
[earth_x,earth_y,earth_z] = sphere(50);
earth_x = earth_x*r_earth + 0;
earth_y = earth_y*r_earth + 0;
earth_z = earth_z*r_earth + 0;
surf (earth_x,-earth_y,earth_z);
h = findobj('Type', 'surface');
%image = 'http://upload.wikimedia.org/wikipedia/commons/thumb/c/cd/Land_ocean_ice_2048.jpg/1024px-Land_ocean_ice_2048.jpg';
image = imread('Earth_Texture.jpg');
imshow(image);
set(h,'CData',image, 'FaceColor','texturemap','edgecolor','none')
grid off;
axis on;
axis equal;
xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');
view([0,0,-1]); % adjust the viewing angle
%zoom(2);
hold off;

r_distance = r_vector_ISS - r_vector_chaser;
r_distance_mod = sqrt(r_distance(1,:).^2+r_distance(2,:).^2+r_distance(3,:).^2);

figure (3);
plot (t/T_chaser, r_distance_mod,'-b');
xlabel ('Nº of Orbits [N]');
ylabel ('Distance ISS-Chaser [km]');


%% C) Calculation of Delta V - Rendezvous
fprintf ("\n\nB) Calculation of Delta V variation for Rendezvous\n");

N_rev_vector = [2:1:30];
delta_V_vector = zeros (length(N_rev_vector),1);
a_chaser_vector = zeros (length(N_rev_vector),1);
ra_chaser_vector = zeros (length(N_rev_vector),1);
N_rev_vector = [2:1:30];

for i= 1:length(N_rev_vector)
    a_chaser = ( (1 - (delta_theta/(2*pi*N_rev_vector(i))) )^(2/3) ) * r_ISS
    e_chaser = (r_ISS/a_chaser) - 1.0 ;
    rp_chaser = a_chaser*(1-e_chaser);
    ra_chaser = a_chaser*(1+e_chaser);
    a_chaser_vector(i) = a_chaser;
    ra_chaser_vector(i) = ra_chaser;

    % Calculation of the Delta V
    va_chaser = sqrt ( ((2*mu)/ra_chaser)-(mu/a_chaser) );
    delta_V = v_ISS - va_chaser;
    delta_V_vector(i) = delta_V;
    
end

% Plot the Graphic
figure (4);
plot (N_rev_vector, delta_V_vector,'-b','LineWidth',4);
xlabel ('Nº of Orbits [N]');
ylabel ('\Delta V [m/s]');

figure (5);
plot (N_rev_vector, a_chaser_vector,'-g','LineWidth',4);
xlabel ('Nº of Orbits [N]');
ylabel ('\Delta V [m/s]');


%% D) Calculation of Delta V - Re-entry
fprintf ("\n\nB) Calculation of Delta V variation for Re-entry\n");

rp_reentry_vector = [r_earth+60:5:r_earth+210];
delta_V_vector_reentry = zeros (length(rp_reentry_vector),1);

for i= 1:length(rp_reentry_vector)
    ra_reentry = r_ISS;
    a_reentry = (ra_reentry+rp_reentry_vector(i))/2;

    % Calculation of the Delta V
    va_reentry = sqrt ( ((2*mu)/ra_reentry)-(mu/a_reentry) );
    delta_V = -va_reentry+v_ISS;
    delta_V_vector_reentry(i) = delta_V;
    
end

% Plot the Graphic
figure (6);
plot (rp_reentry_vector-r_earth, delta_V_vector_reentry,'-b','LineWidth',4);
title('\Delta V Required to Lower Perigee Altitude')
xlabel ('Perigee Altitude of Re-entry Orbit [km]');
ylabel ('(-) \Delta V [km/s]   (opposite direction to velocity)');












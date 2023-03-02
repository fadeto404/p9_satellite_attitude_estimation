clc;
clear;
close all;

%% Costants an parameters
G = 6.674e-11;
EarthRad = 6371*1000;       % m
EarthMass = 5.97e24;        % kg
SatAltitude = 431*1000;     % m
r = EarthRad + SatAltitude;
SatVel = sqrt(G*EarthMass/(r));
SatAngVel = SatVel/r;
SatOrbitTime = 2*pi/SatAngVel;
% xy is the equator then this is the angle of the orbit wrt the equator.
theta = 30;                
%% Draw Earth
[x,y,z] = sphere(50);
surf(x*EarthRad,y*EarthRad,z*EarthRad);
axis equal;
%% Sat Orbit
t = (0:SatAngVel/1000:2*pi);
x = r*cos(t);
y = r*sin(t)*cos(deg2rad(theta));
z = -r*sin(t)*sin(deg2rad(theta));
hold on;
scatter3(x, y, z, Marker=".");
%% IGRFMAGM

%% Animation
simTime = SatOrbitTime*20;
Ts = 0.1/2;
t = (0:Ts:simTime*SatAngVel);
x_initial = r*cos(0);
y_initial = r*sin(0)*cos(deg2rad(theta));
z_initial = -r*sin(0)*sin(deg2rad(theta));
a = scatter3(x_initial, y_initial, z_initial,'MarkerEdgeColor','r','LineWidth',1,Marker="o");
[lat,long,rad] = cart2sph(x_initial,z_initial,y_initial);
[XYZ,H,D,I,F,DXDYDZ,DH,DD,DI,DF] = igrfmagm(rad,lat,long,decyear(datetime('today')));
MagNorm = XYZ/norm(XYZ);
b = plot3([0,MagNorm(1)*1e6],[0,MagNorm(2)*1e6],[0,MagNorm(3)*1e6],'color','k','linewidth',2);
for i = t
    x = r*cos(i);
    y = r*sin(i)*cos(deg2rad(theta));
    z = -r*sin(i)*sin(deg2rad(theta)); 
    drawnow;
    rho = norm([x y z]);
    thetaE = acos(z/rho);
    phiE = 0;
    psiE = atan2(y,x);
    lat = 90-thetaE*180/pi;
    long = psiE*180/pi;
    alt = (rho - EarthRad);
%     [lat,long,rad] = cart2sph(x,y,z);
    [XYZ,H,D,I,F,DXDYDZ,DH,DD,DI,DF] = igrfmagm(alt,lat,long,decyear(datetime('today')));
    BNED = [XYZ(1);XYZ(2);-XYZ(3)];
    BI = TIB(phiE,thetaE+pi,psiE)*BNED;
    MagNorm = BI/norm(BI);
    delete(a);
    delete(b);
    a = scatter3(x, y, z,'MarkerEdgeColor','r','LineWidth',1,Marker="o");
    b = plot3([x,x+MagNorm(1)*1e6],[y,y+MagNorm(2)*1e6],[z,z+MagNorm(3)*1e6],'color','k','linewidth',2);
end
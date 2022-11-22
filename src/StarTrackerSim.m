%% Star Tracker Simulation
clear; clc; close all;
%  rng(1);                         % Consstant RNG
%% Generate star catalog: Uniform random points on a sphere
% Muller method
num_stars = 5000;
u = randn([num_stars, 1]);
v = randn([num_stars, 1]);
w = randn([num_stars, 1]);

unit = symunit;
r = 25*unit.ly; % Distance to stars: 25 lightyears
[r,~] = separateUnits(unitConvert(r,unit.m));
r = double(r);
% r = 10000;

P_norm = sqrt(u.*u + v.*v + w.*w);
P = r*[u./P_norm, v./P_norm, w./P_norm];

figure(1); hold on;
title('Star Map')
scatter3(P(:,1), P(:,2), P(:,3), Marker=".")
grid on; grid minor;
axis equal;
%% True Attitude
% The Optical axis is the z axis of the satelite frame.
% The generated quaternion is the rotation from the global frame to the
% local satellite frame.
x_true = compact(randrot)
x_eul = rad2deg(quat2eul(x_true));
camera_body_axis = [0 0 1];
optic_axis = quatrotate(x_true,camera_body_axis);
camera_mid = optic_axis*r;
x_new = quatrotate(x_true,[1 0 0]);
y_new = quatrotate(x_true,[0 1 0]);

% Plot The Optic Axis
x = [0 , camera_mid(1)];
y = [0 , camera_mid(2)];
z = [0 , camera_mid(3)];
ox = plot3(x,y,z,'color','r','linewidth',2,'DisplayName','Optical Axis');
legend(ox,  'Location','northwest');
%% Star Tracker Specifications
fov = 15;           % Field of View   
focal_l = 1;        % in meters
ppi = 300;
ppm = ppi*100/2.54;     % pixels per m
resolution = 1024;
sensSize = resolution/ppm;
focal_l = sensSize/(2*tan(deg2rad(fov/2)));
KCent = 0.1;
Err = sensSize*KCent/resolution;

rad_fov = r* tan(deg2rad(fov/2));        % Radius of the circle in the fov and r meters away.
photo_scale = focal_l/r;
sqr_len = rad_fov*sqrt(2);
rad_lens = focal_l* tan(deg2rad(fov/2));
lens_edge = rad_lens*sqrt(2);

%Find the edges of the spherical quadrilateral
corners = [lens_edge lens_edge;lens_edge -lens_edge;-lens_edge -lens_edge;-lens_edge lens_edge];
corners = corners./2;
for i = 1:4
    CrnrVec(i,:) = [corners(i,1) corners(i,2) focal_l];
    CrnrVec(i,:) = CrnrVec(i,:)/norm(CrnrVec(i,:));
%     LocalToCrnr_rot(:,:,i) = RU([0 0 1],CrnrVec(i,:));
%     LocalToCrnr(i,:) = rotm2quat(LocalToCrnr_rot(:,:,i));
    LocalToCrnr_ax(i,:) = vrrotvec([0 0 1],CrnrVec(i,:));
    LocalToCrnr(i,:) = axang2quat(LocalToCrnr_ax(i,:));
    GlobalToCrnr(i,:) = quatmultiply(LocalToCrnr(i,:),x_true);
    CrnrAxis(i,:) = quatrotate(GlobalToCrnr(i,:),[0 0 1]);
    CrnrPoints(i,:) = CrnrAxis(i,:)*r;
    CrnrAxis_Lcl(i,:) = quatrotate(LocalToCrnr(i,:),[0 0 1]);
end
% plot3([0,CrnrPoints(1,1)],[0,CrnrPoints(1,2)],[0,CrnrPoints(1,3)],'color','k','linewidth',2);
% plot3([0,CrnrPoints(2,1)],[0,CrnrPoints(2,2)],[0,CrnrPoints(2,3)],'color','g','linewidth',2);
% plot3([0,CrnrPoints(3,1)],[0,CrnrPoints(3,2)],[0,CrnrPoints(3,3)],'color','g','linewidth',2);
% plot3([0,CrnrPoints(4,1)],[0,CrnrPoints(4,2)],[0,CrnrPoints(4,3)],'color','g','linewidth',2);


%% Star Tracker Snap
%Find all the stars insie the FOV of the tracker.
count = 0;
for i = 1:num_stars
    vec = P(i,:)/norm(P(i,:));
    angle =  acos(vec(1)*optic_axis(1) + vec(2)*optic_axis(2) + vec(3)*optic_axis(3));
    if(angle < deg2rad(fov/2))
        count = count + 1; 
        d = dot(camera_mid,optic_axis)/dot(optic_axis,P(i,:));
        %Project the stars on a plane
        projected(count,:) = P(i,:)*d;
        StarsInFOV(count,:) = P(i,:);
        % Move the projected points to the center
        move_to_zero(count,:) = projected(count,:) - camera_mid;
    end
end
for i = 1:4
    d = dot(camera_mid,optic_axis)/dot(optic_axis,CrnrPoints(i,:));
    %Project the stars on a plane
    projectedCrnr(i,:) = CrnrPoints(i,:)*d;
    % Move the projected points to the center
    moveCrnrToZero(i,:) = projectedCrnr(i,:) - camera_mid;
end
scatter3(projected(:,1), projected(:,2), projected(:,3), 'MarkerEdgeColor','k','DisplayName','Projected Stars',Marker=".");
% legend('Projected Points')
scatter3(move_to_zero(:,1), move_to_zero(:,2), move_to_zero(:,3), 'MarkerEdgeColor',"#A2142F",'DisplayName','Centered Stars',Marker=".");
% legend('Centered Points')
% scatter3(moveCrnrToZero(:,1), moveCrnrToZero(:,2), moveCrnrToZero(:,3), 'MarkerEdgeColor','r',Marker=".");

%Rotate the projecte plane into xy axis(MAKE IT 2D)
u = cross(optic_axis,[ 0 0 1]);
ang = atan2(norm(cross(optic_axis,[0 0 1])),dot(optic_axis,[0 0 1]));
rotquat = quaternion(axang2quat([u (ang)]));
In2DPlane = rotatepoint(rotquat,move_to_zero);
CrnrIN2D = rotatepoint(rotquat,moveCrnrToZero);

figure(2); hold on;

% scatter(CrnrIN2D(:,1), CrnrIN2D(:,2), CrnrIN2D(:,3), 'MarkerEdgeColor','r',Marker=".");
title('Stars in 2D')
%Find the angle for boresight rotation
x_point = rotatepoint(rotquat,x_new);
plot([0 x_point(1)*rad_fov*0.8],[0 x_point(2)*rad_fov*0.8],'b--');
scatter(In2DPlane(:,1), In2DPlane(:,2), 'MarkerEdgeColor','k',Marker=".");
legend('Local Frame X Axis','Stars')
rot_angle = atan2(-x_point(2),x_point(1));
Rot_mat = [cos(rot_angle) -sin(rot_angle); sin(rot_angle) cos(rot_angle)];

% Get The Final Snapshot.
n = 0;
for i = 1:count
%     if(In2DPlane(i,1)<sqr_len/2 && In2DPlane(i,1)>-sqr_len/2 && In2DPlane(i,2)<sqr_len/2 && In2DPlane(i,2)>-sqr_len/2)
        n = n+1;
        snapshot(n,:) = In2DPlane(i,1:2);
        % Rotate the Points with the boresight rotation.
        snapshot_rot(n,:) = Rot_mat*snapshot(n,:)';
        % Scale the and flip the photo.(Pinhole camera Model)
        snap(n,:) = snapshot_rot(n,:)* -photo_scale;
        StarsInFrame(n,:) = StarsInFOV(i,:);
%     end
end
for i = 1:4
    % Rotate the Points with the boresight rotation.
    crnr_rot(i,:) = Rot_mat*CrnrIN2D(i,1:2)';
    % Scale the and flip the photo.(Pinhole camera Model)
    crnr(i,:) = crnr_rot(i,:)* -photo_scale;
end
figure(3); hold on;
title('Final SnapShot');
scatter(snap(:,1), snap(:,2), 'MarkerEdgeColor','k',Marker=".");
%% Estimate the Star Tracker Attitude.

for i = 1:n
    % Get The Unit Vector From The Local ody Frame To the Star.
    FocalPlaneToStar(i,:) = [snap(i,1) snap(i,2) -focal_l];
    FocalPlaneToStar(i,:) = -FocalPlaneToStar(i,:)/norm(FocalPlaneToStar(i,:));
%     % Get The Quaternion Rotation from The Local reference To The Star.    
%     LocalToStar_ax(i,:) = vrrotvec([0 0 1],FocalPlaneToStar(i,:));
%     LocalToStar(i,:) = axang2quat(LocalToStar_ax(i,:));
%     
%     GlobalToStar(i,:) = quatmultiply(LocalToStar(i,:),x_true);
%     StarAxis(i,:) = quatrotate(GlobalToStar(i,:),[0 0 1]);
%     StarPoint(i,:) = StarAxis(i,:)*r;
%     % Get The Quaternion Rotation from The Global reference To The Star.  
    StarsInFrameAxis(i,:) = StarsInFrame(i,:)/norm(StarsInFrame(i,:));
%     GlobalToStar_ax(i,:) = vrrotvec([0 0 1],StarsInFrameAxis(i,:));
%     GlobalToStar(i,:) = axang2quat([GlobalToStar_ax(i,1:3) -GlobalToStar_ax(i,4)]);
%     StarPoint2(i,:) = quatrotate(GlobalToStar(i,:),[0 0 1]);
%     StarPoint2(i,:) = StarPoint2(i,:)*r;
%     % Estimate Attitude
%     GlobalToLocal(i,:) = quatmultiply(GlobalToStar(i,:),quatinv(LocalToStar(i,:))); 
% %     GlobalToLocal(i,:) = quatmultiply(quatinv(LocalToStar(i,:)),GlobalToStar(i,:)); 
%     X_est_eul(i,:) = rad2deg(quat2eul(GlobalToLocal(i,:)));
%     Optic_axis_new(i,:) = quatrotate(GlobalToLocal(i,:),camera_body_axis);
end
R_true = quat2rotm(x_true);
Eul_true = quat2eul(x_true,'XYZ');
for i = 1:floor(n/3)
    StarsLocal = [FocalPlaneToStar(3*i-2,:)' FocalPlaneToStar(3*i-1,:)' FocalPlaneToStar(3*i,:)'];
    StarsGlobal = [StarsInFrameAxis(3*i-2,:)' StarsInFrameAxis(3*i-1,:)' StarsInFrameAxis(3*i,:)'];
    R = (StarsLocal/StarsGlobal);
    QuatFin(i,:) = rotm2quat(R);
%     EulCalc(i,:) = quat2eul(QuatFin(i,:),'XYZ');
%     EulErr(i,:) =  EulCalc(i,:) - Eul_true;
end
floor(n/3)
QAvg = avg_quaternion_markley(QuatFin)'
EulErr =  quat2eul(QAvg,'XYZ') - Eul_true

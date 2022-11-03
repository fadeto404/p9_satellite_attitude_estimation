classdef StarTrackerV2
    %Star Tracker With estimated measurmnets from a Star Catalog
    properties
        KCent           {mustBeNumeric}       % Error in pixels
        StarMap         {mustBeNumeric}       % Star Catalog
        SizeStarMap     {mustBeNumeric}       % Num of Stars in Star Catalog
        Dist            {mustBeNumeric}       % Distance to the stars
        Fov             {mustBeNumeric}
        Ppi             {mustBeNumeric}
        Res             {mustBeNumeric}
    end
    
    methods
        function obj = StarTrackerV2(KCent)
            obj.KCent = KCent;
            obj = obj.GenerateStarMap();
            obj.Fov = 15;                     % Focal Length
            obj.Ppi = 300;                    % Camera Sensor PPI
            obj.Res = 1024;                   % Resolution of Sensor.
        end
        
        function obj = GenerateStarMap(obj)
            % Muller method
            obj.SizeStarMap = 10000;
            u = randn([obj.SizeStarMap, 1]);
            v = randn([obj.SizeStarMap, 1]);
            w = randn([obj.SizeStarMap, 1]);
            unit = symunit;
            r = 25*unit.ly; % Distance to stars: 25 lightyears
            [r,~] = separateUnits(unitConvert(r,unit.m));
            obj.Dist = double(r);
            P_norm = sqrt(u.*u + v.*v + w.*w);
            obj.StarMap = obj.Dist*[u./P_norm, v./P_norm, w./P_norm];
        end
        
        function [q_meas,q_err] = MeasureAttitude(obj,q_true)
            % Calc Values of diff parameters.
            optic_axis = quatrotate(q_true,[0 0 1]);
            camera_mid = optic_axis*obj.Dist;
            ppm = obj.Ppi*100/2.54;
            SensorSize = obj.Res/ppm;
            focal_l = SensorSize/(2*tan(deg2rad(obj.Fov/2)));
            rad_fov = obj.Dist* tan(deg2rad(obj.Fov/2));        % Radius of the circle in the fov and r meters away.
            photo_scale = focal_l/obj.Dist;
            sqr_len = rad_fov*sqrt(2);
            Err = SensorSize*obj.KCent/obj.Res;
            % Find all the Stars in FOV and project them on a plane
            count = 0;
            StarsInFOV = zeros(100,3);
            move_to_zero = zeros(100,3);
            for i = 1:obj.SizeStarMap
                vec = obj.StarMap(i,:)/norm(obj.StarMap(i,:));
                angle =  acos(vec(1)*optic_axis(1) + vec(2)*optic_axis(2) + vec(3)*optic_axis(3));
                if(angle < deg2rad(obj.Fov/2))
                    count = count + 1;
                    d = dot(camera_mid,optic_axis)/dot(optic_axis,obj.StarMap(i,:));
                    %Project the stars on a plane
                    projected = obj.StarMap(i,:)*d;
                    StarsInFOV(count,:) = obj.StarMap(i,:);
                    % Move the projected points to the center
                    move_to_zero(count,:) = projected - camera_mid;
                end
            end
            StarsInFOV = StarsInFOV(1:count,:);
            move_to_zero = move_to_zero(1:count,:);
            %Rotate the projecte plane into xy axis(MAKE IT 2D)
            u = cross(optic_axis,[ 0 0 1]);
            ang = atan2(norm(cross(optic_axis,[0 0 1])),dot(optic_axis,[0 0 1]));
            rotquat = quaternion(axang2quat([u (ang)]));
            In2DPlane = rotatepoint(rotquat,move_to_zero);
            %Find the angle for boresight rotation
            x_new = quatrotate(q_true,[1 0 0]);
            x_point = rotatepoint(rotquat,x_new);
            rot_angle = atan2(-x_point(2),x_point(1));
            Rot_mat = [cos(rot_angle) -sin(rot_angle); sin(rot_angle) cos(rot_angle)];
            % Get The Final Snapshot.
            n = 0;
            LocalToStar = zeros(100,3);
            GlobalToStar = zeros(100,3);
            for i = 1:count
                if(In2DPlane(i,1)<sqr_len/2 && In2DPlane(i,1)>-sqr_len/2 && In2DPlane(i,2)<sqr_len/2 && In2DPlane(i,2)>-sqr_len/2)
                    n = n+1;
                    snapshot = In2DPlane(i,1:2);
                    % Rotate the Points with the boresight rotation.
                    snapshot_rot = Rot_mat*snapshot';
                    % Scale the and flip the photo.(Pinhole camera Model)
                    snap = snapshot_rot* -photo_scale + randn(1,2)*Err;
                    LocalToStar(n,:) = [snap(1) snap(2) -focal_l];
                    LocalToStar(n,:) = -LocalToStar(n,:)/norm(LocalToStar(n,:));
                    GlobalToStar(n,:) = StarsInFOV(i,:)/norm(StarsInFOV(i,:));
                end
            end
            GlobalToStar = GlobalToStar(1:n,:);
            LocalToStar = LocalToStar(1:n,:);
            QuatFin = zeros(50,4);
            for i = 1:floor(n/3)
                StarsLocal = [LocalToStar(3*i-2,:)' LocalToStar(3*i-1,:)' LocalToStar(3*i,:)'];
                StarsGlobal = [GlobalToStar(3*i-2,:)' GlobalToStar(3*i-1,:)' GlobalToStar(3*i,:)'];
                R = (StarsLocal/StarsGlobal);
                QuatFin(i,:) = rotm2quat(R);
            end
            QuatFin = QuatFin(1:floor(n/3),:);
            q_meas = avg_quaternion_markley(QuatFin)';
            q_err =  q_true - q_meas;
        end
    end
end


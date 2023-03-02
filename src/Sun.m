%--------------------------------------------------------------------------
%
% Sun: Computes the Sun's geocentric position using a low precision 
%      analytical series
%
% Input:
%   Mjd_TT    Terrestrial Time (Modified Julian Date)
% 
% Output:     
%   rSun      Solar position vector [m] with respect to the 
%             mean equator and equinox of J2000 (EME2000, ICRF)
%
% Last modified:   2015/08/12   Meysam Mahooti
% 
%--------------------------------------------------------------------------
function rSun = Sun (date)

global const % Astronomical Constants

fid = fopen('eop19620101.txt','r');
eopdata = fscanf(fid,'%i %d %d %i %f %f %f %f %f %f %f %f %i',[13 inf]);
fclose(fid);

MJD_UTC = Mjday(date(1),date(2),date(3),date(4),date(5),date(6));
[x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = IERS(eopdata,MJD_UTC,'l');
[UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC] = timediff(UT1_UTC,TAI_UTC);
Mjd_TT = MJD_UTC + TT_UTC/86400;

ep  = const.Rad*(84381.412/3600);     % Obliquity of J2000 ecliptic
T   = (Mjd_TT-const.MJD_J2000)/36525; % Julian cent. since J2000

% Mean anomaly, ecliptic longitude and radius
M = const.pi2 * (0.9931267 + 99.9973583*T)-floor((0.9931267 + 99.9973583*T));                    % [rad]
L = const.pi2 * (0.7859444 + M/const.pi2 + ...
                  (6892*sin(M) + 72.0*sin(2*M))/1296e3) - floor(0.7859444 + M/const.pi2 + ...
                  (6892*sin(M) + 72.0*sin(2*M))/1296e3); % [rad]
r = 149.619e9 - 2.499e9*cos(M) - 0.021e9*cos(2*M);             % [m]


C = cos(-ep);
S = sin(-ep);
R_x = zeros(3,3);

R_x(1,1) = 1.0;  R_x(1,2) =    0.0;  R_x(1,3) = 0.0;
R_x(2,1) = 0.0;  R_x(2,2) =      C;  R_x(2,3) =   S;
R_x(3,1) = 0.0;  R_x(3,2) = -1.0*S;  R_x(3,3) =   C;

% Equatorial position vector
rSun = R_x * [r*cos(L), r*sin(L), 0]'; 


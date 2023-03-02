function MEAS=starmeas(qt,FOV,MTH,sigma,nmax)
%
% function MEAS=starmeas(qt,FOV,MTH,sigma)
%
% This program is used to create a star tracker measurements with 3 sigma noise 
% For certain FOV size (may not be SQUARE) and a certain Magnitude Threshould
% ******************************************************************************
%
% The inputs to this program are:
% 1- The attitude quaternion (m x 3) ; where m is the number of time steps
% 2- The sensor Field of View				(deg)
% 3- MTH = Magnitude THreshold				(min=-1.5, max=6.4)
% 4- The standard deviation or the sensor precision	(rad); default= 17e-6 rad.
% *******************************************************************************
%
% The output of this program will be the body and the intertial star measurements
% The size of the measeurement matrix MEAS is (10*nmax X m), where 
% MEAS = [BT BM IT Av] ; BT, BM are the true and meas. body vectors
% nmax= maximum number of allowable measurement stars 
%
% This program is written by *** Malak A. Samaan *** on Oct. 6 2000.
% *******************************************************************************

cw=180/pi;

nmax=nmax;									% maximum number of allowable measurement stars

tetax=FOV/cw;							% x FOV side (rad)
tetay=FOV/cw;							% y FOV side (rad)

ctet=cos(tetay/2); 
stet=sin(tetay/2);
c2eps=(cos(tetax)+cos(tetay))/(cos(tetay)+1);
ceps=sqrt(c2eps); seps=sqrt(1-c2eps);

s1s=[-stet; +ctet*seps;ctet*ceps ];
s2s=[+stet; +ctet*seps;ctet*ceps ];
s3s=[+stet; -ctet*seps;ctet*ceps ];
s4s=[-stet; -ctet*seps;ctet*ceps ];

%********* Read Star Map ********************************************************

load mappar;
IND=find(MAGN<=MTH);						% Indices of the observable star map
MSR=VS(:,IND);								% Extracts the observable star map

CMIN=s1s'*s3s;								% Minimum interstar cosine (ctmax=s1s'*s3s=s2s'*s4s;) 
cfov2=cos(acos(CMIN)/2);				    % Minimum interstar with respect the Optical Axis

% Sensor Attitude (Body Frame = CS * Sensor Frame)

CS=eye(3);									% Now unity matrix, but then must be established.
OAS=CS(:,3);								% Sensor Optical Axis

% Evaluation of the sensor side-corners (Body Frame)

s1=CS*s1s; s2=CS*s2s; s3=CS*s3s; s4=CS*s4s; 

% Evaluation of the normals of the sensor sides (Body Frame)

stx=sin(tetax); n23=cross(s2,s3)/stx; n41=cross(s4,s1)/stx;
sty=sin(tetay); n12=cross(s1,s2)/sty; n34=cross(s3,s4)/sty;


[m,nnn]=size(qt);						% The number of time steps
MEAS=zeros(m,100);
j500=0;

for j=1:m

% Display when every 500th point is reached
if (j==m), 
 disp(sprintf('ST: Done'))
 j500=0;
end
% j500=j500+1;

   A=quat_att_mat(qt(j,:));						% Calculate The direction cosine matrix 
   n=0;   S=[];
   VI=[];BI=[];
   PS=(OAS'*A)*MSR;
   INDG=find(PS>=cfov2);
   if isempty(INDG)==0,
      for i=INDG,
         vi=MSR(:,i);					% stars in the Inertial frame
         bi=A*vi;							% stars in the Body frame
 %       sip=bi+randn(3,1)*sigma;sip=sip/norm(sip);	
 %       stars in the body frame with measeurment error
 %       Shuster Noise Model, JAS Vol. 38, No. 3, 1990, Eq.(10)
         z1=(bi(1)/bi(3));
         z2=(bi(2)/bi(3));
         rz=sigma^2/(1+z1^2+z2^2)*[(1+z1^2)^2 (z1*z2)^2;(z1*z2)^2 (1+z2^2)^2];
         [u_z,e_z]=eig(rz);
         
         noise1=sqrt(e_z(1,1))*randn(1);
         noise2=sqrt(e_z(2,2))*randn(1);
         noise_corr=u_z*[noise1;noise2];

         alpm=z1+noise_corr(1);
         betam=z2+noise_corr(2);
         sip=[alpm;betam;1]/norm([alpm;betam;1]);
         
         if (n12'*sip<0 & n23'*sip<0 & n34'*sip<0 & n41'*sip<0 & n <nmax),
            n=n+1; 
            S(:,n)=sip;
            VI(:,n)=vi; 
            BI(:,n)=bi;
         end
      end
   end
   av=[ones(1,n) zeros(1,nmax-n)];
   BT=[BI(:)' zeros(1,3*(nmax-n))];			% True Body Measurement
   BM=[S(:)' zeros(1,3*(nmax-n))];			% Actual Body Measerment
   IT=[VI(:)' zeros(1,3*(nmax-n))];  		% True Inertial Measurement
   
   MEAS(j,:)=[BT BM IT av];
%   MEAS(j,:)=IT;
   
end

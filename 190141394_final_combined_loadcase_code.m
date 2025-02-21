clear; clc;

i = 1;
j = 3;
k = 9;
l = 4;

% Point A shall be at the pin support at the left end of the beam. Point B 
% shall be at the roller support at the other end of the beam
 
% 'Av' represents vertical reaction force due to the support at A. The same
% shall go for 'Bv'
  
% The load acting upon load case 1 shall be called 'load1' and its
% horizontal distance from A shall be called 'loaddist1'. The same shall go
% for 'load3'

% 'Span' shall represent the total length of the beam, the horizontal
% distance from point A to B

% 'UDL' shall represent the uniformly distributed load acting on the beam
% in load case 2. 'UDLspan' shall represent the area of beam over which the
% unifromly distributed load acts. It is by deffinition half of the span of
% the beam
   
load1 = 10;
load3 = 25;
loaddist1 = 1+(0.1 * j);
loaddist3 = 1.4 - (0.1 * l);
span = 3+(0.2*i);
UDL = (1 + k);
UDLspan = (span / 2);

% No horizontal force resolution shall be necessary throughout the 
% calculations as all forces act perpendicular to the beam

% When finding moments in order to calculate Av and Bv, moments shall
% always be taken clockwise about A

% The point about which bending moments shall be found shall be called
% point P (though this is irrelevant to calculations) and the bending
% moment about P shall be called 'y'

% The horizontal distance between P and A shall be called 'x'

% Sagging shall be treated as positive in terms of bending moments, as 
% previously stipulated in our MOLE test taken for the 3rd of February

                                                       % Load Case 1
                                                       
% Finding reaction force Av and Bv:
                                                       
Bv = (load1 * loaddist1) / span;
   
% Resolve vertically

Av = load1 - Bv;

% Finding bending moments

x1 = 0:0.1:loaddist1;
 y1 = (Av * x1); 

x2 = loaddist1:0.1:span;
 y2 = (Av * x2) - (load1 * (x2 - loaddist1));
 
x3 = [x1,x2];
y3 = [y1,y2];
 
subplot(2,2,1)
plot(x3,y3,'g')
title ('Load Case 1')
 
xlabel ('Distance from A (m)')
ylabel ('Bending moment (Nm)')
 
subplot(2,2,4)
plot(x3,y3,'g')
hold on
  
xlabel ('Distance from A (m)')
ylabel ('Bending moment (Nm)')

                                                       % Load Case 2

% Finding reaction force Av and Bv:

Bv = ((UDL * UDLspan)*(span-(UDLspan/2))) / span; 

% Resolve vertically

Av = (UDL * UDLspan) - Bv; 

% Finding bending moments:

x1 = 0:0.1:UDLspan;
 y1 = (Av * x1); 
   
x2 = UDLspan:0.1:span; 
 y2 = (Av * x2) - (((UDL * (x2-UDLspan)) .* ((x2-UDLspan)/2)));
 
x3 = [x1,x2];
y3 = [y1,y2];
 
subplot(2,2,2)
plot(x3,y3,'b')
title ('Load Case 2')
 
xlabel ('Distance from A (m)')
ylabel ('Bending moment (Nm)')
 
subplot(2,2,4)
plot(x3,y3,'b')
hold on
 
xlabel ('Distance from A (m)')
ylabel ('Bending moment (Nm)')

                                                       % Load Case 3
                                                       
% Finding reaction force Av and Bv:

Bv = ((load3 * loaddist3) - (load3 * (span - loaddist3))) / span;
  
% Resolve vertically
   
Av = - Bv;

% Finding bending moments:

x1 = 0:0.1:loaddist3;
 y1 = (Av * x1); 

x2 = loaddist3:0.1:(span - loaddist3);
 y2 = (Av * x2) - (load3 * (x2 - loaddist3));
 
x3 = (span - loaddist3):0.1:span;
 y3 = (Av * x3) - (load3 * (x3 - loaddist3)) + (load3 * (x3 - (span - loaddist3)));
 
x4 = [x1,x2,x3];
y4 = [y1,y2,y3];
 
subplot(2,2,3)
plot(x4,y4,'r')
title ('Load Case 3')
 
xlabel ('Distance from A (m)')
ylabel ('Bending moment (Nm)')
 
subplot(2,2,4)
plot(x4,y4,'r')
hold on
 
xlabel ('Distance from A (m)')
ylabel ('Bending moment (Nm)')

                                                       % Load Case 4
                                                       
% Finding reaction force Av and Bv:

Bv = ( ( (UDL * UDLspan) * (UDLspan + UDLspan/2) ) - (load3 * (span - loaddist3) ) + (load1 * loaddist1) + (load3 * loaddist3) ) / span;
   
  % Resolve vertically

Av = load1 + (UDL .* UDLspan) - Bv;

% Finding bending moments:
  
if loaddist1 < UDLspan

     if loaddist1 > loaddist3
    
x1 = 0:0.1:loaddist3;
 y1 = (Av * x1); 
    
x2 = loaddist3:0.1:loaddist1;
 y2 = (Av * x2) - (load3 * (x2 - loaddist3));
 
x3 = loaddist1:0.1:UDLspan;
 y3 = (Av * x3) - (load3 * (x3 - loaddist3)) - (load1 * (x3 - loaddist1));
 
x4 = UDLspan:0.1:(span - loaddist3);
 y4 = (Av * x4) - (load3 * (x4 - loaddist3)) - (load1 * (x4 - loaddist1)) - ((UDL * (x4 - UDLspan)) .* ((x4 - UDLspan)/2)); 
 
x5 = (span - loaddist3):0.1:span;
 y5 = (Av * x5) - (load3 * (x5 - loaddist3)) - (load1 * (x5 - loaddist1)) - ((UDL * (x5 - UDLspan)) .* ((x5 - UDLspan)/2)) + (load3 * (x5 - (span - loaddist3)));

x6 = [x1,x2,x3,x4,x5];
y6 = [y1,y2,y3,y4,y5];
 
subplot(2,2,4)
plot(x6,y6,'m')
title ('Combined Load Cases')
hold on
 
xlabel ('Distance from A (m)')
ylabel ('Bending moment (Nm)')

      % The above code deals with a situation whereby load3 is closer to A than
      % load1 is
 
     else
         
      % The code below deals with a situation whereby load1 is closer to A than
      % load3 is
   
     x1 = 0:0.1:loaddist1;
      y1 = (Av * x1); 
    
     x2 = loaddist1:0.1:loaddist3;
      y2 = (Av * x2) - (load1 * (x2 - loaddist1));
 
     x3 = loaddist3:0.1:UDLspan;
      y3 = (Av * x3) - (load1 * (x3 - loaddist1)) - (load3 * (x3 - loaddist3));
 
     x4 = UDLspan:0.1:(span - loaddist3);
      y4 = (Av * x4) - (load1 * (x4 - loaddist1)) - (load3 * (x4 - loaddist3)) - ((UDL * (x4 - UDLspan)) .* ((x4 - UDLspan)/2)); 
 
     x5 = (span - loaddist3):0.1:span;
      y5 = (Av * x5) - (load1 * (x5 - loaddist1)) - (load3 * (x5 - loaddist3)) - ((UDL * (x5 - UDLspan)) .* ((x5 - UDLspan)/2)) + (load3 * (x5 - (span - loaddist3)));
     
x6 = [x1,x2,x3,x4,x5];
y6 = [y1,y2,y3,y4,y5];

subplot(2,2,4)
plot(x6,y6,'m')
title ('Combined Load Cases')
hold on
 
xlabel ('Distance from A (m)')
ylabel ('Bending moment (Nm)')

     end

% The code above deals with a situation whereby load1 is on A's half of the
% beam
 
else  
    
% The code below deals with a situation whereby load1 is on B's half of the
% beam

     if loaddist1 < (span-loaddist3)
    
x1 = 0:0.1:loaddist3;
 y1 = (Av * x1); 
    
x2 = loaddist3:0.1:UDLspan;
 y2 = (Av * x2) - (load3 * (x2 - loaddist3));
 
x3 = UDLspan:0.1:loaddist1;
 y3 = (Av * x3) - (load3 * (x3 - loaddist3)) - ((UDL * (x3 - UDLspan)) .* ((x3 - UDLspan)/2));
 
x4 = loaddist1:0.1:(span - loaddist3);
 y4 = (Av * x4) - (load3 * (x4 - loaddist3)) - ((UDL * (x4 - UDLspan)) .* ((x4 - UDLspan)/2)) - (load1 * (x4 - loaddist1)); 
 
x5 = (span - loaddist3):0.1:span;
 y5 = (Av * x5) - (load3 * (x5 - loaddist3)) - ((UDL * (x5 - UDLspan)) .* ((x5 - UDLspan)/2)) - (load1 * (x5 - loaddist1)) + (load3 * (x5 - (span - loaddist3))) ;
 
x6 = [x1,x2,x3,x4,x5];
y6 = [y1,y2,y3,y4,y5];

subplot(2,2,4)
plot(x6,y6,'m')
title ('Combined Load Cases')
hold on
 
xlabel ('Distance from A (m)')
ylabel ('Bending moment (Nm)')

     % The code above deals with a situiation whereby load1 is closer to A than
     % the component of load3 that is closest to B
 
     else
   
     % The code below deals with a situation whereby load1 is closer to B than
     % the component of load3 that is closest to B
    
     x1 = 0:0.1:loaddist3;
      y1 = (Av * x1); 
    
     x2 = loaddist3:0.1:UDLspan;
      y2 = (Av * x2) - (load3 * (x2 - loaddist3));
 
     x3 = UDLspan:0.1:(span - loaddist3);
      y3 = (Av * x3) - (load3 * (x3 - loaddist3)) - ((UDL * (x3 - UDLspan)) .* ((x3 - UDLspan)/2));
 
     x4 = (span - loaddist3):0.1:loaddist1;
      y4 = (Av * x4) - (load3 * (x4 - loaddist3)) - ((UDL * (x4 - UDLspan)) .* ((x4 - UDLspan)/2)) + (load3 * (x4 * (span - loaddist3))); 
 
     x5 = loaddist1:0.1:span;
      y5 = (Av * x5) - (load3 * (x5 - loaddist3)) - ((UDL * (x5 - UDLspan)) .* ((x5 - UDLspan)/2)) + (load3 * (x5 - (span - loaddist3))) - (load1 * (x5 - loaddist1));   

x6 = [x1,x2,x3,x4,x5];
y6 = [y1,y2,y3,y4,y5];

subplot(2,2,4)
plot(x6,y6,'m')
title ('Combined Load Cases')
hold on
 
xlabel ('Distance from A (m)')
ylabel ('Bending moment (Nm)')

     end
end

lgd = legend ('Load Case 1','Load Case 2','Load Case 3','Load Case 4');
lgd.FontSize = 15;

% It's unfortunate that the code for the combined load case 4 is so long
% but it was my desire to have the case work for all values of I,J,K, and L
% I was simply not able to figure out an easier way of doing this

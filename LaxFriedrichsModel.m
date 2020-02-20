clear all
close all 
clc 

global k 
global h
global Cx 
global Cy   
global TotSize
global vm
global pm 

global Ex
global NoOfExits

%----------------Variable-------Assignment---------------------------------
k=0.12; %Time step 
h=0.1; %Size Step
a = 1;
x=-a:h:a; 
y=-a:h:a;

vm=0.4; %Vm Could be changed dependant on the type of person 
%0.75m/s - walking speed 
%1.5m/s - run
pm=1; 

tMax=50; %Maximum Time

Size = size(x,2)-1; %Maximum defined for grid calculations
TotSize = size(x,2); 

mid = ceil(TotSize/2); 


%---------------Grid--Initialisation---------------------------------------

P = zeros(size(x,2),size(x,2),floor(tMax/k)); %P(xpos,ypos,time) in matrix form 

%This is where each exit point is defined

% %Positions Exits in midpoints of each wall
 %E1 = [mid,TotSize]; %Exit Positions
 %E2 = [mid,1];  %REMEMBER TO CHANGE NOOFEXITS

 E1 = [TotSize,mid];
 E2 = [1,mid] ;
 
%You cannot set an exit position in a corner! 

NoOfExits = 2;

% Cx = zeros(size(x,2),size(x,2)); 
% Cy = zeros(size(x,2),size(x,2)); 

Ex = zeros(NoOfExits,2);

Ex(1,:) = E1; %Assigning exits to the exit matrix
Ex(2,:) = E2;
% Ex(3,:) = E3; 
% Ex(4,:) = E4; 
% 
% Ex(5,:) = E5; 
% Ex(6,:) = E6; 
% Ex(7,:) = E7; 
% Ex(8,:) = E8; 
%Assigning Exit Position to the matrix 


[Cx,Cy] = DirectionMatrix(Ex, vm, TotSize); %Generates the Directional Vector-Field
%%



%------CFL--Condition-----------------------------------------------------
if abs((2*vm*k)/h) > 1 
msgbox('Simulation Unstable');

end 

%--------------Initial Conditions------------------------------------------
 
for i=1:size(x,2) 
    for j=1:size(x,2)

        P(i,j,1) = 1 - x(i)^2 - y(j)^2; %Simple 2D paraboloid
        
        if P(i,j,1) < 0 
            P(i,j,1) = 0; 
        end 
        
        if isreal(P(i,j,1)) ~= 1 
            P(i,j,1) = 0; 
        end 
        

    end 
end

soom = 0; 

for i=1:size(x,2) %Performed to calculate Phi
    for j=1:size(x,2) 
        soom = soom + P(i,j,1); 
    end 
end 





%%

Xmid = ceil(TotSize/2);
Ymid = ceil(TotSize/2);

%-------------LaxFriedrichs----------------------------------

% 


for m=1:floor(tMax/k) %600
    disp([int2str(m), ' / ' , int2str(floor(tMax/k))]);       
    for i=2:Size   %i is y and j is x       
        for j=2:Size               
            P(i,j,m+1) = 0.25*(P(i+1,j,m) + P(i-1,j,m) + P(i,j+1,m) + P(i,j-1,m))...
                - k/(2*h)*(alpha(P,i+1,j,m,1) - alpha(P,i-1,j,m,1))...
                - k/(2*h)*(alpha(P,i,j+1,m,0) - alpha(P,i,j-1,m,0));
               %B = 0 <- x 
               %B = 1 <- y             
        
        end
    end

   %Copies Boundary Conditions on each edge
   
   %Top, Bottom
   P(1,2:Size,m+1) = P(1,2:Size,m); 
   P(Size+1,2:Size,m+1) = P(Size+1,2:Size,m); 
   
   %Sides, Left, Right 
   P(2:Size,1,m+1) = P(2:Size,1,m); 
   P(2:Size,Size+1,m+1) = P(2:Size,Size+1,m); 
   
   %Corners of Domain 
   P(1,1,m+1) = 0.5*(P(2,1,m+1) + P(1,2,m+1)); 
   P(1,Size+1,m+1) = 0.5*(P(1,Size,m+1) + P(2,Size+1,m+1)); 
   P(Size+1,1,m+1) = 0.5*(P(Size,1,m+1) + P(Size+1,2,m+1)); 
   P(Size+1,Size+1,m+1) = 0.5*(P(Size+1,Size,m+1) + P(Size,Size+1,m+1)); 
   
   %Exit Positions 
   
   for V=1:NoOfExits
       TpEx = Ex(V,:);
       
       if TpEx(1) == 1
           %Y Boundary Top
           P(TpEx(1),TpEx(2),m+1) = P(TpEx(1)+1,TpEx(2),m)...
               - k/(2*h)*(alpha(P,TpEx(1)+1,TpEx(2),m,1) - 0);
       elseif TpEx(1) == TotSize
           %Y Boundary Bottom
           P(TpEx(1),TpEx(2),m+1) = P(TpEx(1)-1,TpEx(2),m)...
               - k/(2*h)*(0 - alpha(P,TpEx(1)-1,TpEx(2),m,1));
       elseif TpEx(2) == 1
           %X Boundary L Side
           P(TpEx(1),TpEx(2),m+1) = P(TpEx(1),TpEx(2)+1,m)...
               - k/(2*h)*(alpha(P,TpEx(1),TpEx(2)+1,m,0) - 0);
       elseif TpEx(2) == TotSize
           %X Boundary R Side
           P(TpEx(1),TpEx(2),m+1) = P(TpEx(1),TpEx(2)-1,m)...
               - k/(2*h)*(0 - alpha(P,TpEx(1),TpEx(2)-1,m,0));
       end
   end
    
    
end 

%---------------Plotting of P----------------------------------------------
%%

XPlot=1:size(x,2); %XPlot is defined for points of P.
YPlot=1:size(x,2);



for t=1:10:floor(tMax/k)
  %t = 1000; 
  surf(x,y,P(YPlot,XPlot,t)); 
  colormap(parula); 
  colorbar; 
  caxis([0 1]); 
  title([int2str(t*k),', vm = 0.9 ']);
  zlim([0 1.2]); 
  xlabel('x'); 
  ylabel('y'); 
  pause(0.25); 
  
  
end



function r=alpha(P,i,j,t,B)
   %B = 0 <- x 
   %B = 1 <- y
   
   global Cx 
   global Cy 
   global pm 
   
   n = 1; 
   
   %Position is (i,j) matrix notation.
   
   if (B==0) %x diff 
        r = (Cx(i,j)*P(i,j,t)*(1-(P(i,j,t)/pm)))^n;        
   elseif (B==1) %y diff  
        r  = (Cy(i,j)*P(i,j,t)*(1-(P(i,j,t)/pm)))^n;
   end          
end 



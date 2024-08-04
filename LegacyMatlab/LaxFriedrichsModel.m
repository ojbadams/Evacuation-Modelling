clear all
close all 
clc 

global Cx 
global Cy 
global pm


%%
%----------------Variable-------Assignment---------------------------------
k=0.01; %Time step 
h=0.1; %Size Step
a = 1; %Domain Size - Always square

x=-a:h:a; 
y=-a:h:a;

vm=3; %Free Speed, if no people how fast people move
pm=1;  %Maximum number of people in 1 grid space. 1 is normalised value

tMax=10; %Maximum Time of iterations

matSize = size(x,2); %Size of Matrix 
mid = floor(matSize/2); %Half Size of Matrix

%%
%%--------------Grid--Initialisation---------------------------------------
P = zeros(size(x,2),size(x,2),floor(tMax/k)); %P(xpos,ypos,time) in matrix form 

E1 = [matSize,mid];
E2 = [1,mid] ;
 
noOfExits = 2; %Change noOfExits depending on how many exits there are

Ex = zeros(noOfExits,2); %exits - y pos stored in first col, x pos in second col

Ex(1,:) = E1; %Assigning exits to the exit matrix
Ex(2,:) = E2;

[Cx,Cy] = DirectionMatrix(Ex, vm, matSize); %Generates the Directional Vector-Field

%%
%%---CFL Condition---------------------------------------------------------

if abs((2*vm*k)/h) > 1 
msgbox('Simulation Unstable');
disp("Calculated Value is " + int2str(2*vm*k/h))
return
end 

%%
%--------------Initial Conditions------------------------------------------

[X,Y] = meshgrid(x,y); 

P(:,:,1) = 1 -X.^2 - Y.^2; %Calculated an upsideown parabola. Any negative values are set to 0.  
P(P(:,:,1) < 0 ) = 0; 

%%
%---------------Lax Friedrichs---------------------------------------------

rangeVec = 2:matSize-1; %This is used to avoid going into out of range 

for m=1:floor(tMax/k) %600
    disp([int2str(m), ' / ' , int2str(floor(tMax/k))]);  
   
    
    P(rangeVec,rangeVec,m+1) = 0.25.*(P(rangeVec+1,rangeVec,m) + ...
        P(rangeVec-1,rangeVec,m) + P(rangeVec,rangeVec+1,m) + P(rangeVec,rangeVec-1,m)) ...
        - k/(2*h).*(alpha(P,rangeVec+1,rangeVec,m,1) - alpha(P,rangeVec-1,rangeVec,m,1)) ...
        - k/(2*h).*(alpha(P,rangeVec,rangeVec+1,m,0) - alpha(P,rangeVec,rangeVec+1,m,0));

    P(1,rangeVec,m+1) = P(2,rangeVec,m+1);
    P(matSize,rangeVec,m+1) = P(matSize-1,rangeVec,m+1); 
    
    P(rangeVec,1,m+1) = P(rangeVec,2,m+1);
    P(rangeVec,matSize,m+1) = P(rangeVec,matSize-1,m+1);     
    
    P(1,1,m+1) = 0.5*(P(2,1,m+1) + P(1,2,m+1)); 
    P(matSize,1,m+1) = 0.5*(P(matSize-1,1,m+1) + P(matSize,2,m+1)); 
    P(1,matSize,m+1) = 0.5*(P(1,matSize-1,m+1) + P(2,matSize,m+1));
    P(matSize,matSize,m+1) = 0.5*(P(matSize,matSize-1,m+1) + P(matSize-1,matSize,m+1));    
end 

%---------------Plotting of P----------------------------------------------
%%

XPlot=1:size(x,2); %XPlot is defined for points of P.
YPlot=1:size(x,2);

for t=1:10:floor(tMax/k)
  surf(x,y,P(YPlot,XPlot,t)); 
  colormap(parula); 
  colorbar; 
  caxis([0 1]); 
  title([num2str(t*k),', vm =  ', num2str(vm)]);
  zlim([0 1.2]); 
  xlabel('x'); 
  ylabel('y'); 
  pause(0.1);
end

%%

function r=alpha(P,yRan,xRan,t,B)
   %B = 0 <- x 
   %B = 1 <- y
   
   global Cx 
   global Cy 
   global pm 
   
   %Position is (i,j) matrix notation.
   
   if B == 0 
       r = Cx(yRan,xRan).*P(yRan,xRan,t).*(1-P(yRan,xRan,t)./pm);
   elseif B == 1 
       r = Cy(yRan,xRan).*P(yRan,xRan,t).*(1-P(yRan,xRan,t)./pm); 
   end             
end 



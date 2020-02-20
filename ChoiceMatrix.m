function cMat = ChoiceMatrix(exits,matSize)

cMat = zeros(matSize,matSize); %This is the choice matrix, each position will be assigned a value from 0...noOfExits#
%cMat assigns values according to the exits matrix. i.e. 1 assigns to
%exits(1,:), 2 assigns to exits(2,:), etc... 

for i=1:matSize
    for j=1:matSize
        temp = sqrt((exits(:,1)-i).^2 + (exits(:,2)-j).^2); %Stores the distance from each position to each exit  
        [~,ind] = min(temp); %finds the minimum value and its position - all we need is the position
        
        cMat(i,j) = ind; %Each exit is assigned a number from 1 2 ... n
        
        %There used to be an IF Statement here to determine if temp matrix
        %had two equidistant exits. 
    end
end
end
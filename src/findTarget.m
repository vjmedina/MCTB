function indexList = findTarget(searchIn, colorTargets)
%FINDTARGET - Finds a series of color strings inside a list of colors,
%             and returns the position (index) of each given color 
%             within that list.
%             Note: This function requires the list elements to be 
%             arranged in columns, each column being a different element
%
% INDEXLIST = FINDTARGET(SEARCHIN, COLORTARGETS)
%
%Author: Victor Medina Heierle
%Last update: 23-Jan-2017
%
%    
    numColors = size(colorTargets,1);
    indexList=[];
    for i=1:numColors
        tmpList = find(searchIn(1,:)==colorTargets(i,1) & searchIn(2,:)==colorTargets(i,2) & searchIn(3,:)==colorTargets(i,3));
        indexList = [indexList tmpList];
    end
end
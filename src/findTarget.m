% Note: This function requires the list elements to be arranged in columns,
% each column being a different element
function indexList = findTarget(searchIn, colorTargets)
    numColors = size(colorTargets,1);
    indexList=[];
    for i=1:numColors
        tmpList = find(searchIn(1,:)==colorTargets(i,1) & searchIn(2,:)==colorTargets(i,2) & searchIn(3,:)==colorTargets(i,3));
        indexList = [indexList tmpList];
    end
end
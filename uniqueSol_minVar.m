function [minVarSolution] = uniqueSol_minVar(allSolution)
varArray =  cell2mat({allSolution.var});
minVar = find(varArray == min(varArray)); 
minVarSolution = allSolution(minVar(1)); % Need to check if two results are also possible

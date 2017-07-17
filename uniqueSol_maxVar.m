function [maxVarSolution] = uniqueSol_maxVar(TaylorSolution)
varArray =  cell2mat({TaylorSolution.var});
maxVar = find(varArray == max(varArray));
maxVarSolution = TaylorSolution(maxVar(1));
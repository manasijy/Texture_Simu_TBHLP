function [UniqueSolution] = uniqueSol_MinPlasticSpin(TaylorSolution)

load('SlipSystem24.mat','SlipSystem'); 
n_T = length(TaylorSolution);
theta = zeros(1,n_T);
for i=1:1:n_T
        
        SS= SlipSystem(TaylorSolution(i).B);
        shear = TaylorSolution(i).xb;
         spin = zeros(3,3);

            for ii=1:1:5, 
                spin = spin + shear(ii)*SS(ii).q.M;
            end
            theta(i) = sqrt(spin(1,2)^2 + spin(1,3)^2 + spin(2,3)^2);
end

mintheta = min(theta); 
q = find(theta == mintheta);
UniqueSolution = TaylorSolution(q(1));% presently only the first solution is selected.
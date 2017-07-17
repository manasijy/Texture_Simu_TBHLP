function [spin,UniqSolution] = calculate_spin(Solution,criteria,SlipSystem)

% SlipSystem = SS_setFCC_fp_function();
% load('FCC_SSSet','SlipSystem');

        switch criteria 
            case 0  %'minvar'
                UniqSolution = uniqueSol_minVar(Solution);
            case 1  %'maxvar'
                UniqSolution = uniqueSol_maxVar(Solution);
            case 2  %'minplasticspin'
                UniqSolution = uniqueSol_MinPlasticSpin(Solution);
        end
    SS= SlipSystem(UniqSolution.B);
    shear = UniqSolution.xb;
    spin = zeros(3,3);
   

    for ii=1:1:5 
        spin = spin + shear(ii)*SS(ii).q.M;
    
    end
  
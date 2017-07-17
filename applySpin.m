function [gmatrix] = applySpin(spin,DCMtrx)

rotdcmat = spin2mat(spin); 
    dcg1=  DCMtrx*rotdcmat';
    gmatrix = rotation('matrix',dcg1); 
end
    
function [rotmatrix] = spin2mat(spin)

rotangle = sqrt(spin(1,2)^2 + spin(1,3)^2 + spin(2,3)^2);
rotaxis = [spin(2,3)/rotangle, spin(3,1)/rotangle, spin(1,2)/rotangle];
r = [rotaxis rotangle];
rotmatrix = vrrotvec2mat(r);
end
function ori = txtfile2ori(filename,cs,ss)
if ~nargin
    name = input('The euler angle file name without .txt extension \n');
    filename = [name, '.txt']; 
end
g = fopen(filename);   
g_matrix = textscan(g, '%f %f %f'); 
fclose(g);
ori = orientation('Euler',[g_matrix{1,1},g_matrix{1,2},g_matrix{1,3}]*degree,cs,ss);
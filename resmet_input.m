function [path] = resmet_input(name,new_gmatrix)

%%
%This fuction takes the output of main_LP.. function and creates a folder
%which contains resmet compatible input file consisting of euler angles at
%each deformation step
%%

[r, c] = size(new_gmatrix);
foldername = [name date];
if exist(foldername,'dir'), foldername = [foldername  num2str(random('Poisson',1))]; end
    
mkdir(foldername);
path = fullfile(cd,foldername); %cd reads the current path in string format
for j=1:1:c
    filename = ['Def_step' num2str(j-1) '.txt'];
    fpath = fullfile(path,filename);
    fid= fopen(fpath,'w');%,'n','US-ASCII');
    fprintf(fid,'This file is created for TexTools software, Resmat Corp \n');
    fprintf(fid,'ODF file name: solodf.HODF \n');
    fprintf(fid,'\nB\t%s\t0 \n',num2str(r));

    for i=1:1:r,fprintf(fid,'\t%4.2f\t%4.2f\t%4.2f\t%d\n',Euler(new_gmatrix(i,j))/degree,1); end
    fclose(fid);
end
end

  
function exportDOS(DOS, filename)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

E=DOS.E;
N=DOS.N;
A=cumsum(N);

fp=fopen(filename, 'w');
fprintf(fp, '%f %f %f\n', [E(:), N(:), A(:)]');
fclose(fp);


end


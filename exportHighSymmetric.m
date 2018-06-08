function exportHighSymmetric(D, file, varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


N=size(D.E,2);
if nargin==3
    M=varargin{1};
else
    M=5;
end
T = D.t;
kx = D.kx;
ky = D.ky;
E = D.E;
p = (length(T)-1)/3;

fp = fopen(file, 'w');

for t=1:3
    for i=max(N/2-M+1,1):min(N/2+M,N)
        %trisurf(T, kx, ky, DE(:,i), 'FaceColor', 'interp','EdgeColor','None');
        for j=p*(t-1)+1:p*t+1
            fprintf(fp, '%f %f %f %f %d\n', T(j), kx(j), ky(j), E(j,i), i);
        end
        fprintf(fp, '\n');
    end
    fprintf(fp, '\n');
end
fclose(fp);

end
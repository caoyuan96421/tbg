function exportDispersion(D, file, varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


N=size(D.E,2);
if nargin==3
    M=varargin{1};
else
    M=3;
end
T = D.t;
B1 = D.B(:,1);
B2 = D.B(:,2);
B3 = 2*B1 + B2;
kx1 = D.kx;
ky1 = D.ky;
kx2 = 1/2 * kx1 - sqrt(3)/2 * ky1;
ky2 = sqrt(3)/2 * kx1 + 1/2 * ky1;
kx3 = 1/2 * kx2 - sqrt(3)/2 * ky2;
ky3 = sqrt(3)/2 * kx2 + 1/2 * ky2;

n=length(kx1);
T=[T;T+n;T+2*n;T+3*n;T+4*n;T+5*n;T+6*n;T+7*n;T+8*n;T+9*n;T+10*n;T+11*n];
kx=[kx1;kx1;-kx1+B3(1);-kx1+B3(1);kx2+B1(1);kx2+B1(1);-kx2+B1(1);-kx2+B1(1);kx3+B1(1);kx3+B1(1);-kx3+B1(1);-kx3+B1(1)];
ky=[ky1;-ky1;ky1+B3(2);-ky1+B3(2);ky2+B1(2);-ky2-B1(2);ky2+B1(2);-ky2-B1(2);ky3+B1(2);-ky3-B1(2);ky3+B1(2);-ky3-B1(2)];
E=repmat(D.E,[12,1]);

fp = fopen(file, 'w');

for i=max(N/2-M+1,1):min(N/2+M,N)
    %trisurf(T, kx, ky, DE(:,i), 'FaceColor', 'interp','EdgeColor','None');
    for t=1:length(T)
        fprintf(fp, '%f %f %f %d\n', kx(T(t,1)), ky(T(t,1)), E(T(t,1), i), i);
        fprintf(fp, '%f %f %f %d\n\n', kx(T(t,2)), ky(T(t,2)), E(T(t,2), i), i);
        fprintf(fp, '%f %f %f %d\n', kx(T(t,1)), ky(T(t,1)), E(T(t,1), i), i);
        fprintf(fp, '%f %f %f %d\n\n\n', kx(T(t,3)), ky(T(t,3)), E(T(t,3), i), i);
    end
end

fclose(fp);

end
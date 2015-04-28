function plotDispersion(D)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


figure;
axis tight;
daspect([1,1,1/2000]);
hold on
N=size(D.E,2);
M=1;
T = delaunay(D.kx, D.ky);
T1 = [];
% Remove zero area triangles due to numerical error
for i=1:size(T,1)
    if norm(cross([D.kx(T(i,1))-D.kx(T(i,2)),D.ky(T(i,1))-D.ky(T(i,2)),0],[D.kx(T(i,2))-D.kx(T(i,3)),D.ky(T(i,2))-D.ky(T(i,3)),0])) > 1e-7
        T1 = [T1 ; T(i,:)];
    end
end
T=T1;
B1 = D.B(:,1);
B2 = D.B(:,2);
B3 = 2*B1 + B2;
kx1 = D.kx;
ky1 = D.ky;
kx2 = 1/2 * kx1 - sqrt(3)/2 * ky1;
ky2 = sqrt(3)/2 * kx1 + 1/2 * ky1;
kx3 = 1/2 * kx2 - sqrt(3)/2 * ky2;
ky3 = sqrt(3)/2 * kx2 + 1/2 * ky2;
%triplot(T,D.kx,D.ky)
cm = jet(2*M);
cc = 0;
for i=max(N/2-M+1,1):min(N/2+M,N)
    cc = cc + 1;
    %surf(D.kx, D.ky, D.E(:,:,i), 'FaceColor', 'interp', 'EdgeColor','None');
    trisurf(T, kx1, ky1, D.E(:,i), 'FaceColor', 'interp','EdgeColor','None');
    trisurf(T, kx1, -ky1, D.E(:,i), 'FaceColor', 'interp','EdgeColor','None');
    trisurf(T, -kx1+B3(1), ky1+B3(2), D.E(:,i), 'FaceColor', 'interp','EdgeColor','None');
    trisurf(T, -kx1+B3(1), -ky1+B3(2), D.E(:,i), 'FaceColor', 'interp','EdgeColor','None');
    trisurf(T, kx2+B1(1), ky2+B1(2), D.E(:,i), 'FaceColor', 'interp','EdgeColor','None');
    trisurf(T, kx2+B1(1), -ky2-B1(2), D.E(:,i), 'FaceColor', 'interp','EdgeColor','None');
    trisurf(T, -kx2+B1(1), ky2+B1(2), D.E(:,i), 'FaceColor', 'interp','EdgeColor','None');
    trisurf(T, -kx2+B1(1), -ky2-B1(2), D.E(:,i), 'FaceColor', 'interp','EdgeColor','None');
    trisurf(T, kx3+B1(1), ky3+B1(2), D.E(:,i), 'FaceColor', 'interp','EdgeColor','None');
    trisurf(T, kx3+B1(1), -ky3-B1(2), D.E(:,i), 'FaceColor', 'interp','EdgeColor','None');
    trisurf(T, -kx3+B1(1), ky3+B1(2), D.E(:,i), 'FaceColor', 'interp','EdgeColor','None');
    trisurf(T, -kx3+B1(1), -ky3-B1(2), D.E(:,i), 'FaceColor', 'interp','EdgeColor','None');
end
hold off
camlight;lighting gouraud;
view(3)

end
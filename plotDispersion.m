function plotDispersion(D)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


figure;
axis tight;
daspect([1,1,1/300]);
hold on
N=size(D.E,3);
M=10;
for i=max(N/2-M+1,1):min(N/2+M,N)
    surf(D.kx, D.ky, D.E(:,:,i));
end
hold off

end
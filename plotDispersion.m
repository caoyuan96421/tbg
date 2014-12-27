function plotDispersion(D)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


figure;
axis tight;
daspect([1,1,10]);
hold on
for i=1:size(D.E,3)
    surf(D.kx, D.ky, D.E(:,:,i));
end

end
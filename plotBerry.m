function plotBerry(D, varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


N=size(D.E,2);
if nargin==2
    M=varargin{1};
else
    M=3;
end

if nargin==3
    index=varargin{2}+N/2;
else
    index=max(N/2-M+1,1):min(N/2+M,N);
end

T = D.t;
kx1 = D.kx;
ky1 = D.ky;
kx2 = 1/2 * kx1 - sqrt(3)/2 * ky1;
ky2 = sqrt(3)/2 * kx1 + 1/2 * ky1;
kx3 = 1/2 * kx2 - sqrt(3)/2 * ky2;
ky3 = sqrt(3)/2 * kx2 + 1/2 * ky2;
B1 = D.B(:,1);
B2 = D.B(:,2);
B3 = 2*B1 + B2;

n=length(kx1);
% T=[T;T+n;T+2*n;T+3*n;T+4*n;T+5*n;T+6*n;T+7*n;T+8*n;T+9*n;T+10*n;T+11*n];
% kx=[kx1;kx1;-kx1+B3(1);-kx1+B3(1);kx2+B1(1);kx2+B1(1);-kx2+B1(1);-kx2+B1(1);kx3+B1(1);kx3+B1(1);-kx3+B1(1);-kx3+B1(1)];
% ky=[ky1;-ky1;ky1+B3(2);-ky1+B3(2);ky2+B1(2);-ky2-B1(2);ky2+B1(2);-ky2-B1(2);ky3+B1(2);-ky3-B1(2);ky3+B1(2);-ky3-B1(2)];
% DE=repmat(D.E,[12,1]);
% DB=repmat(D.Berry,[12,1]);
DE=D.E;
DB=D.Berry;

figure
axis tight;
daspect([1,1,1/2000]);
hold on


for i=index
    trisurf(T, kx1, ky1, DE(:,i), DB(:,i), 'FaceColor', 'interp','EdgeColor','None');
end

hold off
if(numel(index) > 1)
    %camlight;lighting gouraud;
    view(3)
else
    view(2);
end

end
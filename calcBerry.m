function D=calcBerry(D)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

N=size(D.E,2);
T = D.t;
n=length(D.kx);

B=zeros(n,N);
C=zeros(n,N);
for j=1:N
    for i=1:length(T)
        p1=squeeze(D.P(T(i,1), :, j));
        p2=squeeze(D.P(T(i,2), :, j));
        p3=squeeze(D.P(T(i,3), :, j));
        th21=dot(p2, p1); th21=th21/abs(th21);
        th32=dot(p3, p2); th32=th32/abs(th32);
        th13=dot(p1, p3); th13=th13/abs(th13);
        phase=angle(th21*th32*th13);
        B(T(i,1),j) =  B(T(i,1),j) + phase;
        B(T(i,2),j) =  B(T(i,2),j) + phase;
        B(T(i,3),j) =  B(T(i,3),j) + phase;
        C(T(i,1),j) =  C(T(i,1),j) + 1;
        C(T(i,2),j) =  C(T(i,2),j) + 1;
        C(T(i,3),j) =  C(T(i,3),j) + 1;
    end
    disp(j);
end


D.Berry = B./C;
end
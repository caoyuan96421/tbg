function [ DOS ] = DOS( D, res )
%DOS Summary of this function goes here
%   Detailed explanation goes here

    function A=area(a,b,c)
        A=0.5*norm(cross([a(:)-b(:);0],[c(:)-b(:);0]));
    end
N=size(D.E,2);
M=10;
lowband = max(N/2-M+1,1);
highband = min(N/2+M, N);
Emin = min(min(D.E(:,lowband:highband)));
Emax = max(max(D.E(:,lowband:highband)));
Ebin = zeros(res,1);
deltaE = (Emax - Emin)/(res-1);
T = D.t;

for i=1:size(T,1)
    tri=T(i,:);
    k=[D.kx(tri)';D.ky(tri)'];
    ak = area(k(:,1),k(:,2),k(:,3));
    for band=lowband:highband
        % For each triangle, approximate with linear bilinear and put it in the
        % corresponding bin in energy scale
        e=D.E(tri,band);
        eavg = sum(e)/3;
        index = round((eavg-Emin)/deltaE) + 1;
        Ebin(index) = Ebin(index) + ak;
        %[emin, mini] = min(e);
        %[emax, maxi] = max(e);
        %midi = 6 / (mini * maxi);
        % E = [x, y, 1] * coef
        %coef = e \ [[k(:,1)',1];[k(:,2)',1];[k(:,3)',1]]; 
        %p = [e(midi)-coef(3);k(1,mini)*(k(2,maxi)-k(2,mini)) - k(2,mini)*(k(1,maxi)-k(1,mini))] \ ...
        %    [coef(1:2)';[k(2,maxi)-k(2,mini),k(1,mini)-k(1,maxi)]];
        
        
    end
end
DOS.E = linspace(Emin,Emax,res);
DOS.N = Ebin;
plot(DOS.E, DOS.N);
end


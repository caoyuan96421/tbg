classdef TBG < handle
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        theta = 0;
        m=1;
        n=1;
        a0=1;
        a;
        A;              % Superlattice lattice vector
        p;         % p(:,:,1) and p(:,:,2) are lattice points on layer 1 and 2 respectively
    end
    
    methods
        function tbg = TBG(a0, n, m)
            tbg.a0 = a0;
            tbg.n = n;
            tbg.m = m;
            tbg.theta = acos((n^2+4*m*n+m^2)/(2*(n^2+n*m+m^2)));
            a1 = [1.5; -0.5*sqrt(3)] * a0;
            a2 = [1.5; 0.5*sqrt(3)] * a0;
            R1 = [cos(tbg.theta/2), -sin(tbg.theta/2); sin(tbg.theta/2), cos(tbg.theta/2)];
            R2 = [cos(tbg.theta/2), sin(tbg.theta/2); -sin(tbg.theta/2), cos(tbg.theta/2)];
            A1 = R1 * ( n*a1 + m*a2);
            A2 = R1 * (-m*a1 + (n+m)*a2);
            tbg.A = [A1, A2];
            tbg.a = [a1, a2];
            
            tbg.p(:,:,1) = find_points(tbg, n, m, R1*a1, R1*a2);
            tbg.p(:,:,2) = find_points(tbg, m, n, R2*a1, R2*a2);
            
            assert(size(tbg.p,1) == n^2+m^2+n*m);
            disp('Twisted bilayer graphene initialized');
        end
        
        function points = find_points(tbg, p, q, a1, a2)
            % P and Q stores all possible points in lattice coordinates
            [P, Q] = meshgrid(-q:p, 0:(p+2*q));
            r = p^2 + q^2 + p*q;
            % The linear transformation transforms the permitted area into
            % a unit square
            P1 = ((p+q) * P + q * Q) ./ r;
            Q1 = (   -q * P + p * Q) ./ r;
            % Now rule out those points that are not in the unit square
            % after transformation
            In = (P1 >= 0) & (P1 < 1)...
            & (Q1 >= 0) & (Q1 < 1);
            % Convert the lattice coordinates back to real coordinates
            % using provided base vector
            points = P(In) * a1' + Q(In) * a2';
            disp('Complete');
        end
        
        function plot(tbg)
            N = tbg.n+tbg.m;
            a = tbg.a0;
            th = tbg.theta;
            R1 = [cos(th/2), -sin(th/2); sin(th/2), cos(th/2)];
            R2 = [cos(th/2), sin(th/2); -sin(th/2), cos(th/2)];
            c1 = 1.5 * a;
            c2 = sqrt(3) / 2 * a;
            c3 = 0.5 * a;
            fh = figure;
            axis equal tight;
            hold on;
            for i=-max(tbg.n, tbg.m)-2:max(tbg.n, tbg.m)+2
                for j=-2:2*N
                    % Convert lattice coordinate to space coordinate
                    x = c1 * double(i + j);
                    y = c2 * double(j - i);
                    p1 = [x;y] + [-a; 0];
                    p2 = [x;y] + [-c3; c2];
                    p3 = [x;y] + [c3; c2];
                    p4 = [x;y] + [a; 0];
                    line1 = [R1*p1, R1*p2, R1*p3, R1*p4];
                    line2 = [R2*p1, R2*p2, R2*p3, R2*p4];
                    plot(line1(1,:), line1(2,:), 'k', line2(1,:), line2(2,:), 'k');
                end
            end
            
            t1 = tbg.A(:,1);
            t2 = tbg.A(:,2);          
            plot([0; t1(1)], [0; t1(2)], 'g', 'LineWidth', 1.5);
            plot([0; t2(1)], [0; t2(2)], 'g', 'LineWidth', 1.5);
            plot([t1(1); t1(1) + t2(1)], [t1(2); t1(2) + t2(2)], 'g', 'LineWidth', 1.5);
            plot([t2(1); t1(1) + t2(1)], [t2(2); t1(2) + t2(2)], 'g', 'LineWidth', 1.5);
            hold off
            figure
            axis equal tight
            hold on
            plot([0; t1(1)], [0; t1(2)], 'g', 'LineWidth', 1.5);
            plot([0; t2(1)], [0; t2(2)], 'g', 'LineWidth', 1.5);
            plot([t1(1); t1(1) + t2(1)], [t1(2); t1(2) + t2(2)], 'g', 'LineWidth', 1.5);
            plot([t2(1); t1(1) + t2(1)], [t2(2); t1(2) + t2(2)], 'g', 'LineWidth', 1.5);
            scatter(tbg.p(:,1,1), tbg.p(:,2,1),36, 'red');
            scatter(tbg.p(:,1,2), tbg.p(:,2,2),36, 'blue');
            hold off;
        end
    end
    
end


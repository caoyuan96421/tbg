function plotState(tbg, D, n)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    t1 = tbg.A(:,1);
    t2 = tbg.A(:,2);    
    
    
    % Move the coordinate such that AA region is centered
    
    ai = inv(tbg.A);
    c1 = tbg.c(:,:,1) * ai';
    c2 = tbg.c(:,:,2) * ai';
    
    C1 = mod(c1 + 0.5, 1) - 0.5;
    C2 = mod(c2 + 0.5, 1) - 0.5;
    
    tc(:,:,1) = C1 * tbg.A';
    tc(:,:,2) = C2 * tbg.A';
    
    n = n + size(D.E, 2)/2;
    
    fprintf(['Plotting eigen state for band energy E=', num2str(D.E(1,n)), '\n']);
    
    p = reshape(D.P(1,:, n), [size(D.E,2)/2, 2]);
    
    maxP = max(abs(p(:)));
    maxX = max(tc(:,1,1));
    

    figure
    
    subplot(1,2,1);
    axis tight
    daspect([1,1,maxP / maxX]);
    title 'Amplitude'
    hold on
%     plot([0; t1(1)], [0; t1(2)], 'g', 'LineWidth', 1.5);
%     plot([0; t2(1)], [0; t2(2)], 'g', 'LineWidth', 1.5);
%     plot([t1(1); t1(1) + t2(1)], [t1(2); t1(2) + t2(2)], 'g', 'LineWidth', 1.5);
%     plot([t2(1); t1(1) + t2(1)], [t2(2); t1(2) + t2(2)], 'g', 'LineWidth', 1.5);
%     plot([(t1(1)+t2(1))/3; (t1(1)+t2(1))*2/3], [(t1(2)+t2(2))/3; (t1(2)+t2(2))*2/3], 'b', 'LineWidth', 1.5);
%     plot([t1(1); t2(1)], [t1(2); t2(2)], 'b', 'LineWidth', 1.5);
    scatter3(tc(:,1,1), tc(:,2,1),zeros(tbg.N*2,1),36, 'r.');
    scatter3(tc(:,1,2), tc(:,2,2),zeros(tbg.N*2,1),36, 'b.');
    
    pa = p.*conj(p);
    
    scatter3(tc(:,1,1), tc(:,2,1), pa(:,1), 30, 'ko');
    scatter3(tc(:,1,1), tc(:,2,1), pa(:,1), 30, pa(:,1), 'o', 'filled');
    
    scatter3(tc(:,1,2), tc(:,2,2), -pa(:,2), 60, 'ko');
    scatter3(tc(:,1,2), tc(:,2,2), -pa(:,2), 60, pa(:,2), 'o', 'filled');
    colorbar;
    
%     plot3(
%     tri1 = delaunay(tc(:,1,1), tc(:,2,1));
%     tri2 = delaunay(tc(:,1,2), tc(:,2,2));
    hold off; 
    
    % Phase
    
    subplot(1,2,2);
    view(2);
    axis tight
    daspect([1,1,maxP / maxX]);
    title 'Phase'
    
    ph = (angle(p))/(2/3*pi);
    
    hold on;
    scatter3(tc(:,1,1), tc(:,2,1), pa(:,1), 30, 'ko');
    scatter3(tc(:,1,1), tc(:,2,1), pa(:,1), 30, ph(:,1), 'o', 'filled');
    
    scatter3(tc(:,1,2), tc(:,2,2), -pa(:,2), 60, 'ko');
    scatter3(tc(:,1,2), tc(:,2,2), -pa(:,2), 60, ph(:,2), 'o', 'filled');
    colorbar;
    
    hold off;
    
%     figure;
%     axis tight
%     daspect([1,1,1]);
%     trisurf(tri1, tc(:,1,1), tc(:,2,1), p(:,1));

    
%     figure;
%     axis tight
%     trisurf(tri2, tc(:,1,2), tc(:,2,2), p(:,2));
end


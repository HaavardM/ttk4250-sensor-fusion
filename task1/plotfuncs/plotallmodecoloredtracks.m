function plotallmodecoloredtracks(xest, probhat, fignum)
    % Creates one plot containing subplots for all modes
    % xest      estimated states (dim(state) x K)
    % probhat   estimated mode probs (M x K)
    % fignum    first figure number
    S = size(probhat, 1);
    figure(fignum); clf;
    for s = 1:S
        subplot(1, S, s);
        z = zeros(size(xest(1, :)));
        col = probhat(s, :); % probability of mode s (inverse of the second mode if M=2)
        surface([xest(1, :);xest(1, :)],[xest(2, :);xest(2, :)],[z;z],[col;col],'facecol','no','edgecol','interp','linew',2);
        hold on; grid on;
        axis('equal')
        title(sprintf('Mode colored track s=%d', s));
    end
end

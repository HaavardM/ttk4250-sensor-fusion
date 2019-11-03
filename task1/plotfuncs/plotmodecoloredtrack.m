function plotmodecoloredtrack(xest, probhat, s, fignum)
    %PLOTMODECOLOREDTRACK
    % xest      estimated states (dim(state) x K)
    % probhat   estimated mode probs (M x K)
    % s         mode to base colors on
    figure(fignum); clf; hold on; grid on;
    z = zeros(size(xest(1, :)));
    col = probhat(s, :); % probability of mode s (inverse of the second mode if M=2)
    surface([xest(1, :);xest(1, :)],[xest(2, :);xest(2, :)],[z;z],[col;col],'facecol','no','edgecol','interp','linew',2);
    axis('equal')
    title(sprintf('Mode colored track s=%d', s));
end

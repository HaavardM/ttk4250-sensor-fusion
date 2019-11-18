function plotcoloredtrack(track, color, figtitle, scaling, fignum)
    figure(fignum); clf; hold on; grid on;
    %z = zeros(size(track(1, :)));
    surface([track(1, :);track(1, :)],[track(2, :);track(2, :)],[scaling*color; scaling*color],[log(color);log(color)],'facecol','no','edgecol','interp','linew',2);
    axis('equal')
    title(figtitle);
end

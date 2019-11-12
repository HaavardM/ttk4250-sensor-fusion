function plotcoloredtrack(track, color, fignum)
    figure(fignum); clf; hold on; grid on;
    %z = zeros(size(track(1, :)));
    surface([track(1, :);track(1, :)],[track(2, :);track(2, :)],[5*color; 5*color],[log(color);log(color)],'facecol','no','edgecol','interp','linew',2);
    axis('equal')
    title('NIS colored track (xy projection)');
end

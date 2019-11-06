function out_filename = printplot(fig, filename, directory)
    if nargin > 2
        out_filename = strcat(directory, filename);
    else
        out_filename = strcat('../latex/plots/', filename);
    end
    set(fig,'Units','inches');
    screenposition = get(fig,'Position');
    set(fig,...
        'PaperPosition',[0 0 screenposition(3:4)],...
        'PaperSize',[screenposition(3:4)]);
    print(out_filename,'-dpdf');

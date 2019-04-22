function imagescc(v,model,fig_title, path_save)
    defval('fig_title','')
    defval('path_save','Fig/')
    imagesc(model.x/1000, model.z/1000, v);
    if isfield(model,'caxis')
        caxis(model.caxis)
    end
    axis equal tight;
    c = colorbar; ylabel(c,'km/s');
    xlabel('km');
    ylabel('km');
    title(fig_title)    
    set(gca,'FontSize',40)
    set(gcf,'PaperPosition', [0 0 20 15]);
    print(path_save, '-depsc');
end
%% wrapping a set of J files into a structure
% and creating a movie showing convergence
function create_convergence_movie(opts)


if opts.tracking == 0
    disp('Try tracking convergence (opts.tracking = 1) next time to make a movie')
else
    fprintf('create_convergence_movie.m is gathering all iterations into a movie... \n')
    
    % setting up video
    vid = VideoWriter([opts.figFolder, '/historyHD.avi']);
    vid.FrameRate=10;
    
    ex_dir = pwd;
    cd(opts.histFolder);
    
    load('J','J');
    JH = J;
    
    %% stacking history
    
    for i=1:J.evalnum
        load(['J_', num2str(i)],'J');
        JH = [JH J];
    end
        
    %%
    
    
    open(vid)
    figure(2);
    fFWI = J.fFWI;
    
    for i=1:J.evalnum
        load(['J_', num2str(i)],'J');
        
        subplot(311);
        imagesc(reshape(J.m.^-0.5,size(J.gFWI)));
        title(['Evaluation number ', num2str(i), '. Velocity']);
        drawnow;
        
        subplot(312);
        imagesc(J.gFWI);
        title('Gradient')
        
        subplot(313);
        fFWI = [fFWI J.fFWI];
        semilogy(fFWI);
        title('|D-D_obs|^2');
        
        frame = getframe(gcf);
        writeVideo(vid,frame)
    end
    close(vid);
    
    cd(ex_dir);
    
end
end
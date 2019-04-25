%% wrapping a set of J files into a structure
% and creating a movie showing convergence
function create_convergence_movie(opts)

ex_dir = pwd;

cd(opts.histFolder);

load('J','J');
JH=J;

for i=1:J.evalnum
    load(['J_', num2str(i)],'J');
    JH = [JH J];
end 


%%
vid = VideoWriter('historyHD.avi');
vid.FrameRate=10;
open(vid)
figure(2);
fFWI = J.fFWI;

for i=1:J.evalnum
    load(['J_', num2str(i)]);
    
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
    %pause(0.5)
    JH = [JH J];
end  
close(vid);

cd(ex_dir);

end
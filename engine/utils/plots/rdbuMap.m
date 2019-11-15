% creates colormap blue-white-red
% (c) Vladimir Kazei, 2019
function a = rdbuMap()
a = [103,0,31;178,24,43;202,0,32;214,96,77;239,138,98;244,165,130;253,219,199;247,247,247;209,229,240;146,197,222;103,169,207;67,147,195;5,113,176;33,102,172;5,48,97]/255;

%a = imresize(a, [100, 3]);
a = imresize(flipud(a), [100, 3], 'bilinear');



end

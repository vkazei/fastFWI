% creates colormap blue-white-red
% (c) Vladimir Kazei, 2019
function a = rdbuMap()
a = zeros(2001,3);
a(:,:) = NaN;
a(1,:) = [0 0 1];
a(1001,:) = [0.95 0.95 0.95];
a(2001,:) = [1 0 0];
a = fillmissing(a,'linear',1);
end
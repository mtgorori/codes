%%
% load tofmap
%%%
figure;
imagesc(projection*1e6);
set(gca,'YDir','normal');
xlabel('Mfq')
ylabel('óMfq')
c = colorbar;
axis square
axis tight
c.Label.String = '[Ęs]';
% exportfig("H:\result\2018_02_27-_ART_validation\projection",'png',[300,300])
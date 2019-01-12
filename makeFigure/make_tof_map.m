%%
% load tofmap
%%%
figure;
imagesc(projection*1e6);
set(gca,'YDir','normal');
xlabel('送信素子')
ylabel('受信素子')
c = colorbar;
axis square
axis tight
c.Label.String = '[μs]';
% exportfig("H:\result\2018_02_27-_ART_validation\projection",'png',[300,300])
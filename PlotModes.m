function PlotModes(S,ind)

tt = linspace(S.t(1),S.t(end),200);
for i=1:length(ind)
    c = S.coefs{ind(i)};
    display(S.mkn(ind(i),:));
    cc = ppval(S.coefs_pp{ind(i)},tt);
    
    % use this block to plot and comet in the complex plane
    plot(real(c),imag(c),'o');
	%xlim = max(abs(real(c)));
	%ylim = max(abs(imag(c)));
	%axis([-xlim xlim -ylim ylim]);
    hold on
    plot(0,0,'k.','markersize',20)
    comet(real(c(1:10:end)),imag(c(1:10:end)));
    %comet(real(cc),imag(cc));
    
    % use this for a plot of the magnitude
    %plot(tt,abs(cc));
    hold off
    pause
end
function plot_stas(stas)
%
% Usage: plot_stas(stas)

if length(size(stas)) == 5
	[NC, num_lags, ~, L, num_colors] = size(stas);
else
	[num_lags, ~, L, num_colors] = size(stas);
	NC = 1;
end
if NC > 4
	print('the stupid figure will not be the right size if you include more than 4 stas')
	NC = 4;
end

cstrings = {'Lum', 'L-M', 'S'};
clrs = {'k', 'k--', 'r', 'g--' ,'b', 'b-.'};
figure
for cc = 1:NC
	if NC == 1
		k = squeeze(stas);
	else
		k = squeeze(stas(cc,:,:,:,:));
	end

	m = max(abs(k), [], 'all');
	% find best lag
	[~, best_lag] = max( var(reshape(k, num_lags, L*L*num_colors),[], 2) );
	
	tkerns = zeros(num_lags, 6);

	for clr= 1:3
		ax1 = subplot(NC, 5, 5*(cc-1)+clr);
		imagesc(squeeze(k(best_lag, :, :, clr)), [-m, m])
		xticklabels([])
		yticklabels([])
		title(cstrings(clr));
		if clr == 1
			ylabel(sprintf('cell %d', cc), "FontSize",14, "FontWeight",'bold')
		end
		axis square
		drawnow;
		% Identify max and minimum points
		[~, linIdx] = max(reshape(k(best_lag, :, :, clr),L*L,1 ));
		[y_max, x_max] = ind2sub([L,L], linIdx);
		tkerns(:, (clr-1)*2+1) = squeeze(k(:, y_max, x_max, clr));
		[~, linIdx] = min(reshape(k(best_lag, :, :, clr),L*L,1 ));
		[y_min, x_min] = ind2sub([L,L], linIdx);
		tkerns(:, (clr-1)*2+2) = squeeze(k(:, y_min, x_min, clr));
		hold on
		plot(x_max, y_max, 'k.')
		plot(x_min, y_min, 'kx')
	end
	ax2 = subplot(NC, 5, 5*(cc-1)+(4:5));
	hold on
	for clr = 1:6
		plot(0:num_lags-1, tkerns(:,clr), clrs{clr}, LineWidth=1 )
	end
	plot([0,num_lags-1], [0,0],'k')
	axis([0, num_lags-1, -m*1.1,m*1.1])
	% the rest just makes this dumb plot the same dumb height as the imagesc
	inner1 = tightPosition(ax1);  % [L B W H] of inner box
	inner2 = tightPosition(ax2);
	p2 = ax2.Position;
	scale = inner1(4)/inner2(4);
	p2(4) = p2(4)*scale;
	off2 = inner2(2) - ax2.Position(2);
	p2(2) = inner1(2) - off2;
	ax2.Position = p2;
end
%colormap(gray)

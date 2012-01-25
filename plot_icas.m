function plot_icas ( x, prefix = "result_")
	nIC = min(size(x))
	nSamples = max(size(x))

	for i=1:nIC
		plot((0:(nSamples-1))/128, x(:,i))
		xlim([0 (nSamples-1)/128])
		print(strcat(prefix, num2str(i), ".png"), "-dpng")
	end
end


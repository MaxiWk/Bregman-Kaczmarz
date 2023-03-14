% plot residuals with thin shaded area between min and max, 
% line plot of median and thick shaded area between 25/75 quantiles
function medianplot_array = plot_minmax_median_quantiles(stat, solvers, array_on_x_axis)

      num_solvers = length(solvers);
      medianplot_array = zeros(1, num_solvers);  % for legend
      
      if size(array_on_x_axis, 2) == 1
         array_on_x_axis = reshape(array_on_x_axis, size(array_on_x_axis, 2), ...
                                                    size(array_on_x_axis, 1));
      end
      
      hold on

      for i = 1:num_solvers

          solver = solvers{i};
          
          h = fill([array_on_x_axis  fliplr(array_on_x_axis)], [log10(stat.maxs(:,i)')  fliplr(log10(stat.mins(:,i))')], solver.plot_params.minmaxcolor,'EdgeColor', 'none');
          set(h,'facealpha', .5)
          h = fill([array_on_x_axis  fliplr(array_on_x_axis)], [log10(stat.quant75s(:,i)')  fliplr(log10(stat.quant25s(:,i))')], solver.plot_params.quantcolor,'EdgeColor', 'none');
          set(h,'facealpha', .5)
          medianplot_array(i) = plot( array_on_x_axis, log10(stat.medians(:,i)), solver.plot_params.linecolor,'LineWidth',2,...
                                      'LineStyle', solver.plot_params.stroke);
      end
      
      % use adaptive 10^ notation on xaxis
      ax = get(gca);
      yticks = ax.YTick;
      yticks = build_short_rounded_array(yticks, 5);
      num_yticks = length(yticks);
      newyticks = cell(num_yticks,1);
      for j=1:num_yticks
          newyticks(j)=cellstr(strcat('10^{',num2str(yticks(j))));
          newyticks(j)=cellstr(strcat(newyticks(j),'}'));
      end
      ylim([floor(log10(min(stat.mins(:)))), ceil(log10(max(stat.maxs(:))))]); 
      set(gca,'ytick',log10(10.^yticks));
      set(gca,'yticklabel',newyticks);
      
      leg = legend('figure()');
      set(leg,'visible','off')  
          
      axis square 
      
      %hold off
      
      
        % requires that arr is an increasing array and max_length > 1
        function rounded_arr = build_short_rounded_array(arr, max_length)
            rounded_arr = floor(min(arr)); 
                    for k=1:length(arr)     
                        if arr(k) > rounded_arr(end)
                            rounded_arr = [rounded_arr, ceil(arr(k))];
                        end
                    end
            while length(rounded_arr) > max_length
                rounded_arr = rounded_arr(1:2:end);
            end
        end
 
end
    
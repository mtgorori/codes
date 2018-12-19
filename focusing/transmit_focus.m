function [focused_rfdata] = transmit_focus(rfdata,target_element,num_sample,delay_time_all,num_echo_receiver)
% focused_rfdata = zeros(num_sample,num_echo_receiver);
rfdata2cat = rfdata(:,1:num_echo_receiver,:);

for jj = 1:length(target_element)
    % ’x‰„ˆ—
    rfdata2cat(num_sample-delay_time_all(jj)+1:end,:,target_element(jj)) = NaN;
    rfdata2cat(isnan(rfdata2cat)) = 0;
    focused_rfdata = cat(1,zeros(delay_time_all(jj),num_echo_receiver,num_echo_receiver),rfdata2cat);
    %          delay_time = delay_time_all(1,target_element(jj));%[sample]
    %          read_range_rfdata = length(delay_time+1:num_sample);
    %          focused_rfdata(1:read_range_rfdata,:) = focused_rfdata(1:read_range_rfdata,:)...
    %                 +  rfdata(delay_time+1:num_sample,1:num_echo_receiver,target_element(jj));%®‘Š‰ÁZ
end

end
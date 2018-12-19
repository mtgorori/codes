function [focused_rfdata] = transmit_focus(rfdata,target_element,num_sample,delay_time_all,num_echo_receiver)
focused_rfdata = zeros(num_sample,num_echo_receiver);
    for jj = 1:length(target_element)
     % ’x‰„ˆ—
         delay_time = delay_time_all(1,target_element(jj));%[sample]
         read_range_rfdata = length(delay_time+1:num_sample);
         focused_rfdata(1:read_range_rfdata,:) = focused_rfdata(1:read_range_rfdata,:)...
                +  rfdata(delay_time+1:num_sample,1:num_echo_receiver,target_element(jj));%®‘Š‰ÁZ
    end
end
ss;
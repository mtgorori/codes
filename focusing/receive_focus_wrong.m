function [focal_signal_total] = receive_focus_wrong(focused_rfdata,target_element,reference_point)
focal_signal_total = 0;
% hilb_focused_rfdata = hilbert(focused_rfdata);
for jj = 1:length(target_element)
    %受信ビームフォーミング（整相加算のため）
    focal_signal_total = focal_signal_total+ ...
        focused_rfdata(reference_point(1,target_element(1,jj)),target_element(1,jj))/length(target_element);
end
% focal_signal_total  = abs(focal_signal_total);
end
function sendTrigger(trg_type,trg_handle)

if strcmp(trg_type, 'ttl')
    out_ = DaqDOut(trg_handle,0,7); % send trigger
    WaitSecs(0.05);
    out_ = DaqDOut(trg_handle,0,0); % reset trigger
    WaitSecs(0.1);
    
elseif strcmp(trg_type, 'serial')
    trg_data = uint8(1);
    IOPort('Write',trg_handle, trg_data);
    WaitSecs(0.05);
    trg_data = uint8(0);
    IOPort('Write',trg_handle, trg_data);
    WaitSecs(0.1);
elseif strcmp(trg_type, 'labjack')
    
    rest_start_trig = 1;
    if rest_start_trig == 1
        lab_put_code(trg_handle,1);
        rest_start_trig = 0;
    end
    WaitSecs(0.15);
else
    warning('No trigger was sent.')
end


end


function [e] = lab_put_code(L,ecode)

try
    
    L.ljudObj.ePut(L.ljhandle, LabJack.LabJackUD.IO.PUT_DIGITAL_BIT, L.chan{ecode}, 1, 0);
    WaitSecs(.002);
    L.ljudObj.ePut(L.ljhandle, LabJack.LabJackUD.IO.PUT_DIGITAL_BIT, L.chan{ecode}, 0, 0);
    
    e = 0;
catch e
end
end
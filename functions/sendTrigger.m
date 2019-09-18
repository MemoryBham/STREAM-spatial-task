function sendTrigger(trg_type,trg_handle)

if strcmp(trg_type, 'ttl')
    out_ = DaqDOut(trg_handle,0,7); % send trigger
    WaitSecs(0.05);
    out_ = DaqDOut(trg_handle,0,0); % reset trigger
    WaitSecs(0.1);
    
elseif strcmp(trg_type, 'serial')
    trg_data = uint8(1);
    IOPort('Write',trg_handle, trg_data);
    
else
    warning('No trigger was sent.')
end


end
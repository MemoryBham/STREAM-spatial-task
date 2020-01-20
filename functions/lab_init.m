%
% Basic command/response example using MATLAB, .NET and the UD driver.
%
% support@labjack.com
%
try
L.ljasm = NET.addAssembly('LJUDDotNet'); %Make the UD .NET assembly visible in MATLAB
L.ljudObj = LabJack.LabJackUD.LJUD;


    %Used for casting a value to a CHANNEL enum
    L.chanType = LabJack.LabJackUD.CHANNEL.LOCALID.GetType;


    %Open the first found LabJack U3.
    [ljerror, L.ljhandle] = L.ljudObj.OpenLabJack(LabJack.LabJackUD.DEVICE.U3, LabJack.LabJackUD.CONNECTION.USB, '0', true, 0);

    %Start by using the pin_configuration_reset IOType so that all
    %pin assignments are in the factory default condition.
    L.chanObj = System.Enum.ToObject(L.chanType, 0); %channel = 0
    L.ljudObj.ePut(L.ljhandle, LabJack.LabJackUD.IO.PIN_CONFIGURATION_RESET, L.chanObj, 0, 0);
    e = 0;
    
    for i = 1:8
        L.chan{i} = System.Enum.ToObject(L.chanType, i-1); %channel = 0
         L.ljudObj.ePut(L.ljhandle, LabJack.LabJackUD.IO.PUT_DIGITAL_BIT, L.chan{i}, 0, 0);
    end
catch e
end
 e
%ljudObj.Close()
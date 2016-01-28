function CycleAverage = MakeCycleAverage( DataVectorIn, sampling_frequency, frequency, MeanCenterFlag )

%% This routine will average a vector across cycles of the stimulus. The
%% stimulus is assumed to be a sinusoid of phase zero.

% sampling_frequency is the number of measurements per second
% frequency is the stimulus frequency in Hz
% when set, MeanCenterFlag mean-centers the data-vector prior to averaging

% Obtain some basic parameters of the analysis

DataLen=length(DataVectorIn);
CycleTime=(1/frequency)*sampling_frequency; % In units of samples
if (CycleTime < 200)
    CycleTime=CycleTime*2;
end
NumCycles=ceil(DataLen/CycleTime);
PartialCycle=mod(DataLen,CycleTime);

% Fit and remove a linear trend from the data vector if there are three
% stimulus cycles or more in the data.

if (NumCycles>2)
    VectorMean=nanmean(DataVectorIn);
    DataVector = detrend_nan(DataVectorIn);
    if (MeanCenterFlag==0)
        DataVector=DataVector+VectorMean;
    end
    
else
    DataVector=DataVectorIn;
end

for i=1:NumCycles
    if (i==1)
        CycleHolder=squeeze(zeros(CycleTime,1));
        CycleHolder(:)=nan;
        CycleHolder(:)=DataVector(1:CycleTime);
        CycleAverageMatrix=[CycleHolder];
    else
        CycleHolder=squeeze(zeros(CycleTime,1));
        CycleHolder(:)=nan;
        if (i==NumCycles)
            if (PartialCycle==0)
                CycleHolder(:)=DataVector(1+CycleTime*(i-1):end);
            else
                CycleHolder(1:PartialCycle)=DataVector(1+CycleTime*(i-1):end);
            end
        else
            CycleHolder(:)=DataVector(1+CycleTime*(i-1):CycleTime*i);
        end
        CycleAverageMatrix=[CycleAverageMatrix CycleHolder];
    end
end

CycleAverage=nanmean(CycleAverageMatrix,2);

end


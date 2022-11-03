clc;
clear;
close all

St1 = StarTrackerV2(1);

for i = 1:300
    x(i,:) = compact(randrot);
    [q,e] = St1.MeasureAttitude(x(i,:));
    qmeas(i,:) = q;
    EulTrue(i,:) = quat2eul(x(i,:),'XYZ');
    EulMeas(i,:) = quat2eul(qmeas(i,:),'XYZ');
    EulErr(i,:) = EulTrue(i,:) - EulMeas(i,:);
end
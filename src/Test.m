clc;
clear;
close all

St1 = StarTrackerV2(1);
tic
for i = 1:1000
    x(i,:) = compact(randrot);
    EulTrue = quat2eul(x(i,:),'XYZ');
    for j = 1:100
        [q,e] = St1.MeasureAttitude(x(i,:));
        EulMeas(j,:) = quat2eul(q,'XYZ');
        EulErr(j,:) = EulTrue - EulMeas(j,:);
    end
    VarEul(i,:) = var(EulErr);
end
toc
mean(VarEul)

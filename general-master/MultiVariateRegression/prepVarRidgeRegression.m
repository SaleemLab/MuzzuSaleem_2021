function [bestInd, spikeRate, traj, posbins, gain, eye_XPos, eye_YPos, pupilSize, lick, smthBallSpd,smthBallAcc,smthTrajSpd,...
    binnedBallSpd,binnedBallAcc, binnedTrajSpd, binnedEye_XPos, binnedEye_YPos, binnedPupilSize, binnedLick] = prepVarRidgeRegression(Nav, spikeTrain, subsetref, nPosbins, Nquantiles)


bestInd =  1:size(spikeTrain,2);

spikeRate = spikeTrain(:,bestInd)';
traj = Nav.trajPercent;
roomlength = 100;%max(floor(traj+1));
traj = normalise1var(traj, roomlength, [], [0 roomlength]);
posbins = 0:roomlength;

gain = Nav.gain;

eye_XPos = Nav.eyeXpos;
eye_YPos = Nav.eyeYpos;
pupilSize = Nav.pupilSize;

pupilSize = (pupilSize - median(pupilSize(subsetref)))./quantile(pupilSize(subsetref),0.95);
eye_XPos  = (eye_XPos - median(eye_XPos(subsetref)))./quantile(eye_XPos(subsetref),0.95);
eye_YPos  = (eye_YPos - median(eye_YPos(subsetref)))./quantile(eye_YPos(subsetref),0.95);
lick = Nav.lick;

smthBallSpd = Nav.smthBallSpd;
smthBallSpd = (smthBallSpd-nanmedian(smthBallSpd(subsetref)))./quantile(smthBallSpd(subsetref),0.95);

smthTrajSpd = Nav.smthBallSpd.*Nav.gain;
smthTrajSpd = (smthTrajSpd-nanmedian(smthTrajSpd(subsetref)))./quantile(smthTrajSpd(subsetref & ~isnan(smthTrajSpd)),0.95);

smthBallAcc = Nav.smthBallSpd;
smthBallAcc = (smthBallAcc-nanmedian(smthBallAcc(subsetref)))./quantile(smthBallAcc(subsetref),0.95);

[binnedBallSpd, Nqout] = linearBinning(smthBallSpd,Nquantiles,subsetref);
[binnedBallAcc, Nqout] = linearBinning(smthBallAcc,Nquantiles,subsetref);
[binnedTrajSpd, Nqout] = linearBinning(smthTrajSpd,Nquantiles,subsetref);
[binnedLick, Nqout] = linearBinning(lick,Nquantiles,subsetref);
[binnedPupilSize, Nqout] = linearBinning(pupilSize,Nquantiles,subsetref);
[binnedEye_XPos, Nqout] = linearBinning(eye_XPos,Nquantiles,subsetref);
[binnedEye_YPos, Nqout] = linearBinning(eye_YPos,Nquantiles,subsetref);
end
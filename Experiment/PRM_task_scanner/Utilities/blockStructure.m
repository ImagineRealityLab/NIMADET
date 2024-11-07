function [blocks,miniblocks] = blockStructure(nOri,nMB,nRep)

ori = 1:nOri;

imaOri = ori(randperm(nOri)); % shuffle which orientation to start with
imaOri = repmat(imaOri,1,nRep);
blocks = imaOri';

perOri = ori(randperm(nOri)); % shuffle which orientation to start with
perOri = repmat(perOri,1,nRep*nMB);
miniblocks = perOri';
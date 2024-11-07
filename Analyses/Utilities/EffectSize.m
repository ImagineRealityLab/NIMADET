function effect_size = EffectSize(M1,M2,SD1,SD2)
% calculates the effect size as cohen's d for t-tests

SDpooled = sqrt((SD1^2+SD2^2)/2);

effect_size = (M2-M1)/SDpooled;
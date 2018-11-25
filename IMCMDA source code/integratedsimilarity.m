function [sd,sm] = integratedsimilarity(FS,FSP,SS,SSP,kd,km)         
sm = FS.*FSP+km.*(1-FSP);            
sd = SS.*SSP+kd.*(1-SSP);            
end
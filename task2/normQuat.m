function [q] = normQuat(qin)
%NORMQUAT Summary of this function goes here
%   Detailed explanation goes here
q = qin / norm(qin);
end


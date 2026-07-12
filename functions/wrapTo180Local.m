% wrapTo180Local.m — Wrap degree angles into (-180, 180].
% Maps any real angle (deg) onto the principal branch without relying on
% Mapping Toolbox wrapTo180. Used for longitude innovations in the UKF
% and for azimuth residuals in batch orbit determination, so that a
% jump across +/-180 deg does not look like a large measurement error.
% Math: wrapped = mod(angle + 180, 360) - 180.
% Author: Pasquale Marzaioli

function wrapped = wrapTo180Local(angle)
% Map degree-valued angles to the principal interval (-180, 180].
wrapped = mod(angle + 180, 360) - 180;
end

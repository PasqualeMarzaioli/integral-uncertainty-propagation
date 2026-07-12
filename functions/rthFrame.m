% rthFrame.m — Radial-transverse-normal (RTN/RTH) basis from state.
% Builds orthonormal row basis that maps ECI vectors into local orbit frame:
%   r_hat = r/|r|,  h_hat = (r x v)/|r x v|,  t_hat = h_hat x r_hat.
% Author: Pasquale Marzaioli

function frame = rthFrame(state)
% Build row basis vectors that map ECI components into radial-transverse-normal.
radial = state(1:3) / norm(state(1:3));
normal = cross(state(1:3), state(4:6));
normal = normal / norm(normal);
transverse = cross(normal, radial);
frame = [radial.'; transverse.'; normal.'];
end

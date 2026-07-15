% rthFrame.m — Build the local radial-transverse-normal rotation.
% Forms an orthonormal row basis from the supplied inertial position and velocity.
% Author: Pasquale Marzaioli

function frame = rthFrame(state)
% Build row basis vectors that map inertial components into radial-transverse-normal.
radial = state(1:3) / norm(state(1:3));
normal = cross(state(1:3), state(4:6));
normal = normal / norm(normal);
transverse = cross(normal, radial);
frame = [radial.'; transverse.'; normal.'];
end

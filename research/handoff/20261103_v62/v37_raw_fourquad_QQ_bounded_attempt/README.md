# Direct exact QQ raw four-quadric attempt: bounded failure

This distinct direct QQ run reads the literal Fraction-polynomial 24-by-56
source frozen in v36. It never lifts reduced coefficients back to QQ. The source
is pinned to that exact native QQ input. It attempts a QQ standard basis,
cokernel annihilator and U-saturation with the same known U as the field pilot.

The single zero-worker Singular process reaches only the input-ready marker.
It is killed at the sixty-second cap (60.022 seconds; maximum RSS 719508 KiB),
before any standard-basis, annihilator or saturation result. Its process group
is terminated and no solver remains. There is no QQ positive or negative
membership/rank conclusion. No longer retry or larger QQ solver is included.

The complete source script, literal source pin, runner, stdout and resource
receipt are retained. The next separate method is a fixed-column homogeneous
finite-field pilot; no QQ claim is inferred from the previously successful
inhomogeneous field saturation.

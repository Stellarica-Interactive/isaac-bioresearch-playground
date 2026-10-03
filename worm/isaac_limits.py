"""Constraints the Isaac stage enforces, importable without starting Isaac.

``worm/isaac/stage.py`` cannot be imported outside a running simulator -- it
pulls in ``pxr`` -- so a constant that both the Isaac articulation and the
quasi-static solver must agree on cannot live there if the second is to be
testable on its own.

Kept deliberately small. Anything here is a number two body solvers share, and a
disagreement between them would make their results incomparable.
"""

from __future__ import annotations

import math

#: Joint limit in degrees, matching ``stage.DEFAULT_JOINT_LIMIT_DEG``. Real worms
#: bend harder than this; the limit exists to stop the chain folding through
#: itself, and ``pinned`` in the runner's report is the fraction of joints
#: against it -- above zero the mechanics rather than the nervous system are
#: setting the body shape.
JOINT_LIMIT_DEG = 60.0

JOINT_LIMIT_RAD = math.radians(JOINT_LIMIT_DEG)

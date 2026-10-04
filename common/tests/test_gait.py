"""Undulation metrics, tested against signals whose answer is known in advance.

The metrics these replace were wrong for four sections of
``docs/model_assumptions.md`` without any test failing, because every test they
had compared them against another unvalidated number. The tests here construct
sinusoids with a chosen period and a chosen inter-joint phase, so there is a right
answer to miss.
"""

from __future__ import annotations

import numpy as np
import pytest

from common.body.gait import dominant_period_samples, measure, travel_at

HZ = 240.0


def _wave(
    *,
    period_s: float,
    phase_deg: float,
    amplitude_deg: float = 20.0,
    joints: int = 24,
    seconds: float = 20.0,
) -> np.ndarray:
    """A travelling wave with a known period and a known phase per joint."""
    t = np.arange(int(seconds * HZ))[:, None]
    j = np.arange(joints)[None, :]
    w = 2 * np.pi / (period_s * HZ)
    return np.deg2rad(amplitude_deg) * np.sin(w * t - np.deg2rad(phase_deg) * j)


def test_a_known_travelling_wave_reads_back_its_own_phase() -> None:
    for phase in (5.0, 15.0, 23.0, 40.0):
        g = measure(_wave(period_s=2.0, phase_deg=phase))
        assert g.inter_joint_phase_deg == pytest.approx(phase, abs=1.5)
        assert g.is_travelling_wave


def test_a_standing_oscillation_has_no_phase() -> None:
    """Every joint moving together, at large amplitude. The failure mode the
    original metric was written to avoid, and did avoid."""
    g = measure(_wave(period_s=2.0, phase_deg=0.0))
    assert g.amplitude_deg > 5.0
    assert abs(g.inter_joint_phase_deg) < 1.0
    assert not g.is_travelling_wave


def test_direction_is_signed() -> None:
    forward = measure(_wave(period_s=2.0, phase_deg=23.0)).inter_joint_phase_deg
    backward = measure(_wave(period_s=2.0, phase_deg=-23.0)).inter_joint_phase_deg
    assert forward > 0 > backward
    assert forward == pytest.approx(-backward, abs=1.0)


def test_a_static_bend_is_not_an_oscillation() -> None:
    bent = np.tile(np.deg2rad(np.linspace(0.0, 40.0, 24)), (4800, 1))
    g = measure(bent)
    assert g.amplitude_deg == pytest.approx(0.0, abs=1e-9)
    assert not g.is_travelling_wave


def test_amplitude_survives_a_period_longer_than_any_fixed_window() -> None:
    """The bug this module was written for.

    The committed metric measured amplitude over the last 240 samples -- one
    second. A body whose dominant period is ten seconds swings just as far, but a
    one-second window catches a sliver of one excursion and reports almost
    nothing. That is how a moving body came to be recorded as ``amp 0.96deg``.
    """
    slow = _wave(period_s=10.0, phase_deg=20.0, amplitude_deg=20.0, seconds=40.0)
    g = measure(slow)
    assert g.period_s == pytest.approx(10.0, abs=0.5)
    assert g.amplitude_deg > 10.0, "a slow swing is still a swing"

    one_second = slow[-int(HZ) :]
    naive = float(np.degrees((one_second - one_second.mean(0)).std()))
    assert naive < 0.5 * g.amplitude_deg, "the old window would have missed it"


def test_the_fixed_lag_understated_every_real_wave() -> None:
    """Closed form: ``travel(lag) = 2 sin(phase) sin(2 pi lag / T)``.

    At the committed 12-sample lag and a 2 s period that factor is 0.156, so a
    genuine gait could not read above about +0.12 however cleanly it propagated.
    This is the arithmetic that made ``travel +0.12`` look like a weak wave.
    """
    w = _wave(period_s=2.0, phase_deg=23.0)
    fluct = w - w.mean(axis=0, keepdims=True)

    period = 2.0 * HZ
    at_committed_lag = travel_at(fluct, 12)
    at_quarter = travel_at(fluct, int(period / 4))

    expected_ratio = np.sin(2 * np.pi * 12 / period)
    assert at_committed_lag / at_quarter == pytest.approx(expected_ratio, abs=0.02)
    assert at_committed_lag < 0.2 < at_quarter, "sixfold compression"


def test_propagation_is_not_reported_when_no_period_resolves() -> None:
    """Refusing to answer is the point.

    A window holding less than two cycles puts the spectral peak in the lowest
    bin regardless of content, so any phase computed from it describes the window
    length. The metric returns nan rather than a number.
    """
    rng = np.random.default_rng(0)
    drift = np.cumsum(rng.standard_normal((600, 24)) * 1e-4, axis=0)
    g = measure(drift)
    assert not np.isfinite(g.period_s)
    assert g.inter_joint_phase_deg == 0.0
    assert not g.is_travelling_wave


def test_too_short_a_recording_is_refused() -> None:
    # Compared field by field rather than against a whole Gait: nan != nan, so
    # dataclass equality is never true for a refusal.
    for too_small in (np.zeros((4, 24)), np.zeros((4800, 1))):
        g = measure(too_small)
        assert g.amplitude_deg == 0.0
        assert not np.isfinite(g.period_s)
        assert g.inter_joint_phase_deg == 0.0
        assert g.travel == 0.0
        assert not g.is_travelling_wave


def test_a_frozen_body_is_not_given_a_phase() -> None:
    """The failure that prompted MIN_AMPLITUDE_DEG.

    Ablating VD froze the body at ``bend 51.1 deg, extent 0.24, amp 0.00 deg`` and
    the metric reported a phase of **+34.8 deg** -- larger than the scripted
    wave's genuine +23.0 deg, from a body that was not moving. The correlation was
    running on PhysX jitter, which is far above the 1e-12 floor in `travel_at`.
    """
    rng = np.random.default_rng(1)
    # A static coil, plus jitter with a genuine phase gradient, so the test fails
    # if the guard is removed rather than merely passing for lack of structure.
    t = np.arange(4800)[:, None]
    j = np.arange(24)[None, :]
    coil = np.deg2rad(np.linspace(0.0, 51.0, 24))[None, :]
    jitter = 1e-7 * np.sin(2 * np.pi * t / 480.0 - np.deg2rad(30.0) * j)
    frozen = coil + jitter + rng.standard_normal((4800, 24)) * 1e-9

    g = measure(frozen)
    assert g.amplitude_deg < 0.1
    assert g.inter_joint_phase_deg == 0.0, "no oscillation, so no phase to report"
    assert g.travel == 0.0
    assert not g.is_travelling_wave

    # The same shape at a real amplitude still reads its phase, so the guard is a
    # floor on noise and not a blanket refusal.
    real = coil + np.deg2rad(10.0) * np.sin(2 * np.pi * t / 480.0 - np.deg2rad(30.0) * j)
    assert measure(real).inter_joint_phase_deg == pytest.approx(30.0, abs=2.0)


# -- is there an oscillation at all? ----------------------------------------
#
# MIN_CYCLES_IN_WINDOW bounds how long a period may be, which is not the same
# question as whether there is one. With no periodicity the spectrum is roughly
# flat, argmax lands on the lowest bin the length check admits -- period exactly
# half the window -- and the check passes it. Observed in a run that reported
# `period 6.0s` at t=12, `7.0s` at t=14 and `8.0s` at t=16, each half the elapsed
# window, on 1.3 degrees of decaying amplitude, with `travel +0.28` derived from
# it. See MIN_PEAK_PROMINENCE and model_assumptions 5AK.


def _travelling(seconds: float, hz: float, *, waves: float = 1.5, amp_deg: float = 16.0):
    """A travelling wave over a given span, for the prominence tests.

    Named apart from the module's own `_wave`: defining a second helper under
    that name silently replaced the first and broke five unrelated tests with a
    TypeError, which is a neater demonstration of why than any comment.
    """
    t = np.arange(int(seconds * hz)) / hz
    j = np.arange(23)
    phase = 2.0 * np.pi * (0.5 * t[:, None] - waves * j[None, :] / 23)
    return np.radians(amp_deg) * np.sin(phase)


def test_a_real_wave_is_still_resolved() -> None:
    """The guard must not cost the thing the metric exists to measure."""
    gait = measure(_travelling(20.0, 240.0), hz=240.0)
    assert gait.period_s == pytest.approx(2.0, abs=0.05)
    assert gait.travel > 0.5


@pytest.mark.parametrize(
    "name",
    ("noise", "drift", "decay"),
)
def test_a_signal_with_no_rhythm_resolves_no_period(name: str) -> None:
    """Each of these used to return a period of exactly half the window.

    `travel` is computed at a lag of a quarter period, so a fabricated period
    produces a fabricated propagation: these reported up to +0.28 for bodies
    doing nothing rhythmic at all.
    """
    rng = np.random.default_rng(0)
    samples, joints = 240 * 20, 23
    t = np.arange(samples) / 240.0
    signals = {
        "noise": np.radians(1.3) * rng.standard_normal((samples, joints)),
        "drift": np.radians(5.0) * np.linspace(0.0, 1.0, samples)[:, None] * np.ones(joints),
        "decay": np.radians(1.3)
        * np.exp(-t / 8.0)[:, None]
        * rng.standard_normal((samples, joints)),
    }
    gait = measure(signals[name], hz=240.0)
    assert not np.isfinite(gait.period_s), f"{name} resolved a period of {gait.period_s}"
    assert gait.travel == 0.0, f"{name} reported travel {gait.travel}"


def test_a_slow_rhythm_just_inside_the_margin_still_resolves() -> None:
    """The case both guards had to keep: a real but slow oscillation.

    Written first as *two* cycles in the window, which was the floor at the time.
    MIN_CYCLES_IN_WINDOW was then raised to 2.5 precisely because a period of
    exactly half the window is where a spectrum with nothing to offer puts its
    peak, so two cycles is now rejected by design and that version of this test
    contradicted the fix. The guard is still about prominence rather than about
    which bin the peak lands in -- it simply needs a margin above the floor.
    """
    hz = 240.0
    samples = int(20.0 * hz)
    t = np.arange(samples) / hz
    # Four cycles in twenty seconds: period 5 s, comfortably inside the margin
    # and still far slower than the 2 s gait the metric was built for.
    slow = np.radians(10.0) * np.sin(2.0 * np.pi * t / 5.0)[:, None] * np.ones(23)
    period = dominant_period_samples(slow - slow.mean(axis=0, keepdims=True))
    assert np.isfinite(period), "a clean four-cycle oscillation must resolve"
    assert period / hz == pytest.approx(5.0, abs=0.3)


def test_exactly_two_cycles_is_rejected_on_purpose() -> None:
    """The cost of the margin, asserted so it is a decision and not a surprise.

    A period of exactly half the window is indistinguishable from the artifact
    of §5AK.3, so it is refused. A real rhythm loses nothing permanently: its
    period stays put as the window grows, so it moves off the boundary and is
    then accepted, while an artifact tracks the boundary forever.
    """
    hz = 240.0
    samples = int(20.0 * hz)
    t = np.arange(samples) / hz
    boundary = np.radians(10.0) * np.sin(2.0 * np.pi * t / 10.0)[:, None] * np.ones(23)
    assert not np.isfinite(dominant_period_samples(boundary - boundary.mean(axis=0, keepdims=True)))
    # The same rhythm in a longer window clears the margin and resolves.
    longer = np.arange(int(40.0 * hz)) / hz
    same = np.radians(10.0) * np.sin(2.0 * np.pi * longer / 10.0)[:, None] * np.ones(23)
    period = dominant_period_samples(same - same.mean(axis=0, keepdims=True))
    assert np.isfinite(period), "four cycles of the same rhythm must resolve"
    assert period / hz == pytest.approx(10.0, abs=0.5)

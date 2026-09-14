"""Whether a body is undulating, and whether the undulation travels.

Both questions are easy to answer with a number that looks reasonable and means
nothing. The metrics this module replaces did exactly that, in two ways that had
gone unnoticed through §5C to §5Q of ``docs/model_assumptions.md``:

* **Amplitude was measured over a one-second window.** For a body whose dominant
  period is ten seconds -- which is what the connectome-driven model actually
  does -- a one-second window sees a fraction of one excursion and reports
  ``0.96°`` for a body that swings ``14.75°``. Slow movement read as no movement.
* **Propagation was measured at a fixed lag of twelve samples**, 50 ms at 240 Hz.
  For two joints oscillating at period ``T`` with phase difference ``φ``, the
  difference-of-correlations formula evaluates in closed form to

      travel(lag) = 2 sin(φ) sin(2π·lag/T)

  (verified against synthetic sinusoids to three decimals). The reading is
  therefore the real signal multiplied by ``sin(2π·lag/T)`` -- a factor set by an
  arbitrary constant, not by the body. At ``T`` = 2 s that factor is **0.156**, so
  a scripted wave which is a genuine gait by construction, covering 2.39 body
  lengths in twenty seconds, reported ``travel +0.12``. That was read as a weak
  but real wave. It was a strong wave seen through a metric compressed sixfold,
  and the docstring's claim that "+1 is a clean head-to-tail wave" was unreachable
  at that lag for any body.

The fix in both cases is to derive the timescale from the signal instead of
assuming one. ``sin(2π·lag/T)`` is maximised at ``lag = T/4``, where the reading is
exactly ``2 sin(φ)`` -- so measuring there both maximises sensitivity and makes the
number invertible into something with units: the phase difference between adjacent
joints, in degrees, which is what "the wave travels" actually means.

This does not overturn §5O, §5P or §5Q. Re-measured correctly, the connectome
baseline's inter-joint phase is **−0.5°** against the scripted wave's **+23.0°**:
still no propagation. What changes is that every magnitude in those sections was
on a compressed scale, and a weak wave would have been invisible.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

#: Longest period, in samples, that is still treated as an oscillation rather than
#: drift. A "period" as long as the recording is not a measurement of anything --
#: the FFT simply puts its peak in the lowest bin. ASSUMED; it only guards the
#: degenerate case.
MIN_CYCLES_IN_WINDOW = 2.0

#: Amplitude below which no phase is reported, in degrees. A body frozen in a
#: bend still jitters by a fraction of a float, and correlating that jitter
#: returns a confident, structured, meaningless number: ablating VD froze the body
#: at ``bend 51.1 deg, extent 0.24, amp 0.00 deg`` and the phase came back
#: **+34.8 deg**, larger than the scripted wave's genuine +23.0 deg. The guard in
#: :func:`travel_at` did not catch it because PhysX jitter is far above 1e-12.
#:
#: ASSUMED. A tenth of a degree is well under any real undulation here -- the
#: scripted wave is 21 deg and even the quietest connectome run is 3 deg -- and
#: well above numerical noise, but nothing measures it.
MIN_AMPLITUDE_DEG = 0.1


@dataclass(frozen=True)
class Gait:
    """What a run of joint angles says about undulation."""

    amplitude_deg: float
    """Temporal standard deviation of joint angle, over whole cycles."""

    period_s: float
    """Dominant oscillation period. ``nan`` if the window holds too few cycles."""

    inter_joint_phase_deg: float
    """Phase difference between adjacent joints. Positive is head to tail.

    This is the quantity that distinguishes crawling from flexing in place. Zero
    means every joint does the same thing at the same moment, which is a standing
    oscillation however large its amplitude.
    """

    travel: float
    """``2 sin(phase)``, the raw reading at the quarter-period lag.

    Kept because it is what the older metric reported and comparisons against the
    committed numbers in ``docs/model_assumptions.md`` need it -- but note those
    were taken at a fixed 12-sample lag and are smaller by ``sin(2π·12/T)``.
    """

    @property
    def is_travelling_wave(self) -> bool:
        """Deliberately strict, and deliberately not a headline number.

        A body must both oscillate appreciably and carry a real phase gradient.
        The threshold is ASSUMED -- there is no measurement behind 5 degrees -- so
        this is a convenience for reading a table, never evidence on its own.

        This was the only thing standing between a frozen, coiled body and a
        reported phase of +34.8 degrees; :data:`MIN_AMPLITUDE_DEG` now stops the
        number being produced at all, which is the right place for the guard.
        """
        return bool(
            np.isfinite(self.period_s)
            and self.amplitude_deg > 1.0
            and abs(self.inter_joint_phase_deg) > 5.0
        )


def dominant_period_samples(fluct: np.ndarray) -> float:
    """Period of the strongest oscillation in the mid-body joint, in samples.

    The mid-body joint rather than a pooled average: averaging across joints that
    are out of phase cancels precisely the signal being looked for.
    """
    if fluct.shape[0] < 4:
        return float("nan")
    mid = fluct[:, fluct.shape[1] // 2]
    spectrum = np.abs(np.fft.rfft(mid))
    spectrum[0] = 0.0  # the mean is already removed; bin 0 is numerical residue
    peak = int(np.argmax(spectrum))
    if peak == 0:
        return float("nan")
    period = fluct.shape[0] / peak
    if period * MIN_CYCLES_IN_WINDOW > fluct.shape[0]:
        return float("nan")
    return float(period)


def travel_at(fluct: np.ndarray, lag: int) -> float:
    """Difference of correlations between adjacent joints at ``+lag`` and ``-lag``.

    A body flexing in place correlates equally in both directions, so the
    difference is zero for anything symmetric in time; a wave moving head to tail
    correlates more strongly at positive lag. Each joint's temporal mean is
    removed by the caller, so a static bend contributes nothing.
    """
    if lag < 1 or 2 * lag >= fluct.shape[0] or fluct.shape[1] < 2:
        return 0.0
    anterior = fluct[lag:-lag, :-1].ravel()

    def corr(shift: int) -> float:
        posterior = fluct[lag + shift : fluct.shape[0] - lag + shift, 1:].ravel()
        if anterior.std() < 1e-12 or posterior.std() < 1e-12:
            return 0.0
        return float(np.corrcoef(anterior, posterior)[0, 1])

    return corr(+lag) - corr(-lag)


def measure(history: list[np.ndarray] | np.ndarray, *, hz: float = 240.0) -> Gait:
    """Measure undulation over whole cycles, at the lag the body itself sets.

    The window is the entire recording rather than a fixed tail: any fixed window
    is either shorter than the period -- which is how ``0.96°`` was reported for a
    body swinging ``14.75°`` -- or an assumption about a timescale that has not
    been measured yet.
    """
    angles = np.asarray(history, dtype=np.float64)
    if angles.ndim != 2 or angles.shape[0] < 8 or angles.shape[1] < 2:
        return Gait(0.0, float("nan"), 0.0, 0.0)

    fluct = angles - angles.mean(axis=0, keepdims=True)
    period = dominant_period_samples(fluct)
    if not np.isfinite(period):
        # No resolvable oscillation. Amplitude is still meaningful; propagation is
        # not, and returning a number for it would be the mistake this module
        # exists to stop.
        return Gait(float(np.degrees(fluct.std())), float("nan"), 0.0, 0.0)

    # Whole cycles only: a partial cycle biases the standard deviation by however
    # much of the swing happens to fall inside the window.
    cycles = int(fluct.shape[0] // period)
    whole = fluct[-int(cycles * period) :]
    whole = whole - whole.mean(axis=0, keepdims=True)

    amplitude = float(np.degrees(whole.std()))
    if amplitude < MIN_AMPLITUDE_DEG:
        # A phase needs an oscillation to be the phase *of*. Below this the
        # correlation is computed on jitter and returns whatever the jitter
        # happens to align to -- see MIN_AMPLITUDE_DEG.
        return Gait(amplitude, float(period / hz), 0.0, 0.0)

    lag = max(int(round(period / 4)), 1)
    travel = travel_at(whole, lag)
    # travel = 2 sin(phase) at the quarter-period lag, so phase inverts directly.
    # Clipped because noise can push |travel| past 2, where the arcsine is
    # undefined -- that is a saturated reading, not a phase above 90 degrees.
    phase = float(np.degrees(np.arcsin(np.clip(travel / 2.0, -1.0, 1.0))))
    return Gait(
        amplitude_deg=amplitude,
        period_s=float(period / hz),
        inter_joint_phase_deg=phase,
        travel=float(travel),
    )

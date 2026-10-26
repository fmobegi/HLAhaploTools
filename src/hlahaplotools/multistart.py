"""EM from several start values, keeping the best fit.

The EM finds a stationary point of the likelihood near its start. With many
loci or rare haplotypes, the default start (haplotype counts over all phases,
``numberOccurrence``) can sit on a symmetric stationary point that is not the
maximum: when the phases of a person involve haplotypes nobody else carries,
an even start stays even and the EM stops there. Perturbing the start breaks
that symmetry. In the benchmark (scripts/benchmark.py multistart), perturbed
starts came much closer to the best fit than fully random ones.

:func:`estimate_multistart` runs the requested start first, then perturbed
starts (``numberOccurrence`` with seeded noise) with seeds ``seed + 1,
seed + 2, ...``, and returns the estimate with the highest log-likelihood, so
results stay reproducible.
"""

from __future__ import annotations

from hlahaplotools._core import estimate


def estimate_multistart(
    dataset, nomenclature, resolutions=None, *, starts: int = 10, extra_init: str = "perturbation", **options
):
    """Run :func:`hlahaplotools.estimate` from ``starts`` start values; keep the best.

    The first start uses ``options["init"]`` (default ``numberOccurrence``) and
    ``options["seed"]`` (default 0); the others use ``extra_init`` (default
    ``perturbation``) with seeds ``seed + 1`` to ``seed + starts - 1``. Returns
    ``(best, summary)`` where ``summary`` lists ``(init, seed, log_likelihood,
    iterations, converged)`` for every start, in the order run.
    """
    if starts < 1:
        raise ValueError("starts must be at least 1")
    seed = options.pop("seed", 0)
    init = options.pop("init", "numberOccurrence")
    best, summary = None, []
    for k in range(starts):
        this_init, this_seed = (init, seed) if k == 0 else (extra_init, seed + k)
        est = estimate(dataset, nomenclature, resolutions, init=this_init, seed=this_seed, **options)
        summary.append((this_init, this_seed, est.log_likelihood, est.iterations, est.converged))
        # Ties keep the earlier start, so the result does not depend on float noise.
        if best is None or est.log_likelihood > best.log_likelihood + 1e-9 * abs(best.log_likelihood):
            best = est
    return best, summary

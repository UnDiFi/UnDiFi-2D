# tests/tools

`compare.py` -- tolerance-aware regression comparator. See its module
docstring for the full format description and rationale (ROADMAP.md
Part II, Phase 0 / issue #11); short version:

```
pip install -r requirements.txt

python3 compare.py scalar <reference> <candidate> [--rtol R] [--atol A] [--json PATH]
python3 compare.py shock  <reference> <candidate> [--rtol R] [--atol A] [--json PATH]
python3 compare.py selftest
```

`../checking_fitting_runs.sh` and `../checking_capturing_runs.sh` call
`compare.py scalar` against the checked-in `checksum_*` files instead of
doing an exact `diff`, so a build/compiler/parallelization change that
only perturbs the last few digits of a converged run no longer registers
as a failure.

`compare.py selftest` exercises both the `scalar` and `shock` comparison
paths against fixtures already in the repo and asserts: an identical copy
passes, a 1-ulp perturbation fails at zero tolerance but passes at the
default tolerance, and a gross perturbation always fails. Run it after
touching this file, or as a quick check that the tolerances still make
sense.

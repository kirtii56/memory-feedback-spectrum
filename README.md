# Memory feedback and the linearized spectrum

Verification code for the note *When Does Memory Feedback Change the Linearized Spectrum?*

## Files
- `meas_verify_final.py` recomputes the numerical and symbolic values reported in Section 2 and prints a pass/fail table (25 checks).
- `check_hyp33.py` constructs explicit trajectory pairs for Hypothesis 3.3' in the five examples (Section 2.6).
- `requirements.txt` needs only SymPy 1.14.0. Everything else is the Python standard library.

## Run
    pip install -r requirements.txt
    python3 meas_verify_final.py
    python3 check_hyp33.py

Expected: `TOTAL 25 checks; 0 FAIL` and `ALL PASS`.

## What this does and does not show
The scripts reproduce illustrative numbers and check algebraic identities. They do not prove the theorems; the proofs are in the paper.

## Runs recorded
| Run | Environment | Result |
|---|---|---|
| 1 | Colab, Python 3.13.15, SymPy 1.14.0 | four documentation values corrected, then re-checked |
| 2 | same image | identical to run 1 |
| 3 | Linux, Python 3.13.16, SymPy 1.14.0 | 25 checks, 0 failures |

## License
Add a LICENSE file of your choice (for example MIT) before publishing.

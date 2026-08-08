#!/usr/bin/env python3
"""Benchmark legacy vs PETSc linear solver backends on the radial-flux PM machine,
sweeping several PETSc configurations and verifying the results agree.

Usage:
    benchmark_radial_machine.py [model...] [options]

Options:
    --iterations N         timed repetitions per backend/config (default 3;
                           results are reported as the median of N)
    --pin CORE             pin the solver to a CPU core (recommended on loaded
                           machines)
    --petsc-options "label:opts"
                           run PETSc with the given PETSC_OPTIONS, labelled for
                           the output table.  May be given multiple times.
                           "label:" with empty opts tests the built-in default.
                           If the argument has no 'label:' prefix, the options
                           themselves are used as the label.
    --tolerance REL        verification tolerance: the PETSc result must agree
                           with legacy to within REL * max|A| for every node
                           (default 1e-3)
    --no-verify            skip the result comparison

The mesh is generated once per model with fmesher and reused for every backend
so the comparison isolates the solver.  Each backend/config is timed on the
same mesh; runs are interleaved to cancel CPU/thermal drift.

Output columns:
    model    config    nodes    legacy (s)    petsc (s)    speedup
             nonlin    maxAbsErr    relErr    check
where relErr is maxAbsErr / max|A| over the field and 'check' is PASS/FAIL
against the requested tolerance (only meaningful when the node ordering of the
two .ans files matches, which it does because the same mesh is used).
"""
import argparse
import glob
import os
import shutil
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))  # repository root
BIN = os.path.join(ROOT, "cfemm", "bin")
DATA = os.path.join(ROOT, "mfemm", "testing", "radial_machine", "data")
WORK = os.path.join(tempfile.gettempdir(), "rmbench_work")

FSOLVER = os.path.join(BIN, "fsolver")
FMESHER = os.path.join(BIN, "fmesher")


def run(cmd, env=None, cwd=None, pin=None):
    if pin is not None:
        cmd = ["taskset", "-c", str(pin)] + list(cmd)
    e = dict(os.environ)
    if env:
        e.update(env)
    return subprocess.run(cmd, env=e, cwd=cwd, capture_output=True, text=True)


def wall(fn):
    t0 = time.perf_counter()
    rc = fn()
    return time.perf_counter() - t0, rc


def parse_configs(arguments):
    """Turn the repeated --petsc-options values into [(label, options), ...]."""
    if not arguments:
        return [("default", "")]
    configs = []
    for value in arguments:
        label, sep, opts = value.partition(":")
        if not sep:
            opts = value
            label = value if value else "default"
        configs.append((label, opts))
    return configs


def parse_solution(data):
    """Parse the [Solution] node rows of an .ans file into (x, y, re, im)."""
    text = data.decode("utf-8", errors="replace")
    lines = [ln.strip() for ln in text.splitlines()]
    try:
        idx = lines.index("[Solution]")
    except ValueError:
        return None
    n = int(lines[idx + 1])
    nodes = []
    for i in range(n):
        parts = lines[idx + 2 + i].split()
        if len(parts) < 4:
            return None
        x = float(parts[0])
        y = float(parts[1])
        re = float(parts[2])
        im = float(parts[3]) if len(parts) >= 5 else 0.0
        nodes.append((x, y, re, im))
    return nodes


def compare_solutions(a_bytes, b_bytes, tol):
    """Return (maxAbsErr, relErr, status) for two solution files."""
    if a_bytes is None or b_bytes is None:
        return None, None, "no-output"
    a = parse_solution(a_bytes)
    b = parse_solution(b_bytes)
    if a is None or b is None:
        return None, None, "parse-failed"
    if len(a) != len(b):
        return None, None, f"node-count {len(a)} vs {len(b)}"
    # the same mesh is reused, so the node ordering must match
    scale = 0.0
    max_abs = 0.0
    for (xa, ya, rea, ima), (xb, yb, reb, imb) in zip(a, b):
        if abs(xa - xb) > 1e-9 or abs(ya - yb) > 1e-9:
            return None, None, "mesh-order"
        va = (rea * rea + ima * ima) ** 0.5
        vb = (reb * reb + imb * imb) ** 0.5
        scale = max(scale, va, vb)
        err = ((rea - reb) ** 2 + (ima - imb) ** 2) ** 0.5
        max_abs = max(max_abs, err)
    rel = max_abs / scale if scale > 1e-12 else max_abs
    status = "PASS" if max_abs <= tol * scale else "FAIL"
    return max_abs, rel, status


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("models", nargs="*")
    ap.add_argument("--iterations", type=int, default=3)
    ap.add_argument("--pin", type=int, default=None,
                    help="CPU core to pin the solver to")
    ap.add_argument("--petsc-options", action="append", default=[],
                    metavar="label:opts",
                    help="PETSC_OPTIONS variant to benchmark (repeatable)")
    ap.add_argument("--tolerance", type=float, default=1e-3,
                    help="verification tolerance relative to max|A|")
    ap.add_argument("--no-verify", action="store_true")
    args = ap.parse_args()

    configs = parse_configs(args.petsc_options)

    if args.models:
        models = args.models
    else:
        models = [os.path.join(DATA, "radial_machine_sliding.fem")]
        models += sorted(glob.glob(os.path.join(DATA, "radial_machine_redraw_*.fem")))

    os.makedirs(WORK, exist_ok=True)

    legacy_env = {"XFEMM_SOLVER_BACKEND": "legacy"}
    petsc_envs = []
    for label, opts in configs:
        env = {"XFEMM_SOLVER_BACKEND": "petsc"}
        if opts:
            env["PETSC_OPTIONS"] = opts
        petsc_envs.append((label, opts, env))

    header = (f"{'model':<32s}{'config':<10s}{'nodes':>7s}"
              f"{'legacy(s)':>11s}{'petsc(s)':>11s}{'speedup':>9s}"
              f"{'nonlin':>7s}")
    if not args.no_verify:
        header += f"{'maxAbsErr':>12s}{'relErr':>10s}{'check':>7s}"
    print(header)
    print("-" * len(header))

    totals = {label: 0.0 for label, _, _ in petsc_envs}
    legacy_total = 0.0

    for fem in models:
        if not os.path.isabs(fem):
            cand = os.path.join(os.getcwd(), fem)
            if not os.path.exists(cand):
                cand = os.path.join(DATA, fem)
            fem = cand
        base = os.path.splitext(os.path.basename(fem))[0]
        meshdir = os.path.join(WORK, "mesh_" + base)
        os.makedirs(meshdir, exist_ok=True)
        # mesh a local copy (fmesher writes mesh files next to the input file)
        local_fem = os.path.join(WORK, os.path.basename(fem))
        with open(fem, "rb") as fin, open(local_fem, "wb") as fout:
            fout.write(fin.read())
        r = run([FMESHER, os.path.basename(fem)], cwd=WORK)
        if r.returncode != 0:
            print(f"{base}: meshing failed\n{r.stdout}\n{r.stderr}")
            continue
        for suffix in (".node", ".ele", ".pbc", ".edge"):
            src = os.path.join(WORK, base + suffix)
            if os.path.exists(src):
                os.replace(src, os.path.join(meshdir, base + suffix))
        with open(os.path.join(meshdir, base + ".node")) as f:
            nodes = f.readline().split()[0]

        def copy_mesh():
            for suffix in (".node", ".ele", ".pbc", ".edge"):
                s = os.path.join(meshdir, base + suffix)
                if os.path.exists(s):
                    d = os.path.join(WORK, base + suffix)
                    if os.path.exists(d):
                        os.remove(d)
                    shutil.copy2(s, d)

        def solve_once(env):
            copy_mesh()
            dt, rc = wall(lambda: run([FSOLVER, base], env=env, cwd=WORK,
                                      pin=args.pin))
            ans = None
            if rc.returncode == 0:
                p = os.path.join(WORK, base + ".ans")
                if os.path.exists(p):
                    with open(p, "rb") as f:
                        ans = f.read()
            return dt, rc, ans

        # Warm up every backend/config once (Cuthill-McKee ordering and, for
        # PETSc, the first factorization fall outside the timed loop).
        solve_once(legacy_env)
        for _, _, env in petsc_envs:
            solve_once(env)

        ltimes = []
        lrc = 0
        lnit = 0
        legacy_ans = None
        ctimes = {label: [] for label, _, _ in petsc_envs}
        crc = {label: 0 for label, _, _ in petsc_envs}
        cnit = {label: 0 for label, _, _ in petsc_envs}
        cans = {label: None for label, _, _ in petsc_envs}

        # Interleave legacy and every PETSc config to cancel CPU/thermal drift.
        for _ in range(args.iterations):
            dt, rc, ans = solve_once(legacy_env)
            ltimes.append(dt)
            lrc |= rc.returncode
            lnit = max(lnit, rc.stdout.count("Newton Iteration"),
                       rc.stdout.count("Successive Approx"))
            if ans is not None:
                legacy_ans = ans
            for label, _, env in petsc_envs:
                dt, rc, ans = solve_once(env)
                ctimes[label].append(dt)
                crc[label] |= rc.returncode
                cnit[label] = max(cnit[label], rc.stdout.count("Newton Iteration"),
                                  rc.stdout.count("Successive Approx"))
                if ans is not None:
                    cans[label] = ans

        ltimes.sort()
        lt = ltimes[len(ltimes) // 2]
        legacy_total += lt

        for label, opts, _ in petsc_envs:
            ctimes[label].sort()
            pt = ctimes[label][len(ctimes[label]) // 2]
            totals[label] += pt
            nit = max(lnit, cnit[label])
            speedup = lt / pt if pt > 0 else 0.0
            row = (f"{base:<32s}{label:<10s}{nodes:>7s}{lt:>11.2f}{pt:>11.2f}"
                   f"{speedup:>8.2f}x{nit:>7d}")
            if not args.no_verify:
                max_abs, rel, status = compare_solutions(
                    legacy_ans, cans[label], args.tolerance)
                if max_abs is None:
                    row += f"{'':>12s}{'':>10s}{status:>7s}"
                else:
                    row += (f"{max_abs:>12.3e}{rel:>10.3e}{status:>7s}")
            if crc[label] or lrc:
                row += f"  [rc legacy={lrc}, petsc={crc[label]}]"
            print(row)
            sys.stdout.flush()

    if len(models) > 1:
        print("-" * len(header))
        for label, opts, _ in petsc_envs:
            speedup = legacy_total / totals[label] if totals[label] > 0 else 0
            print(f"{'TOTAL':<32s}{label:<10s}{'':>7s}{legacy_total:>11.2f}"
                  f"{totals[label]:>11.2f}{speedup:>8.2f}x")


if __name__ == "__main__":
    main()

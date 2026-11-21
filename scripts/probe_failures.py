#!/usr/bin/env python3
"""Explore the region of parameter space where MESA fails to build an
initial (pre-MS) model, or to evolve it to the ZAMS, in EM.

Each parameter set is run in its own process (bin/em_probe), because a
failure inside create_star is a hard abort. Parameters are handled in
the AMP scaled units used by ff_amp2 (each in [0,1]):

    M     = d1 + 0.75
    Z     = 10**(0.8*d2 - 2.1)
    Y     = 0.1*d3 + 0.22
    alpha = 2*d4 + 1

Results are appended to a CSV cache (--out), so interrupted or repeated
explorations reuse earlier runs.

Modes:
    point   run the base point only
    scan    1-D scans through the base point along each axis
    edges   bisect along each axis (both directions) for the
            pass/fail boundary nearest the base point
    grid    2-D grid through the base point in two chosen axes
    tc      try each of a list of pre-MS central temperatures (--tc-values)
            on its own, at each point of the scan neighbourhood or of a
            sweep along one axis (--axis, --range); prints a point x T_c table

Examples:
    probe_failures.py point --scaled 0.31 0.52 0.88 0.47
    probe_failures.py scan  --scaled 0.31 0.52 0.88 0.47 --width 0.05 -n 11 -j 8
    probe_failures.py edges --physical 1.06 0.308 0.0089 1.94 --tol 1e-3
    probe_failures.py grid  --scaled 0.31 0.52 0.88 0.47 --axes M Y -n 9 -j 8
"""

import argparse
import concurrent.futures as cf
import csv
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
import threading
import time

AXES = ['M', 'Z', 'Y', 'alpha']

# Outcome codes and the symbols used in text maps

OK, CREATE_FAIL, ZAMS_FAIL, LOGG_FAIL, CRASH, TIMEOUT = \
    'ok', 'create_fail', 'zams_fail', 'logg_fail', 'crash', 'timeout'
SYMBOL = {OK: '.', CREATE_FAIL: 'C', ZAMS_FAIL: 'Z', LOGG_FAIL: 'G', CRASH: '!', TIMEOUT: 'T'}

FIELDS = ['stage', 'd_M', 'd_Z', 'd_Y', 'd_alpha', 'M', 'Y', 'Z', 'alpha', 'f_ov', 'max_age', 'T_c',
          'outcome', 'T_c_used', 't_code', 'exit_code', 'runtime',
          'zams_age', 'zams_Teff', 'zams_L', 'zams_R', 'zams_Xc', 'zams_model',
          'logg_age', 'logg_Teff', 'logg_L', 'logg_R', 'logg_Xc', 'logg_model', 'logg_t_code', 'message']

ZAMS_FIELDS = ['zams_age', 'zams_Teff', 'zams_L', 'zams_R', 'zams_Xc', 'zams_model']
LOGG_FIELDS = ['logg_age', 'logg_Teff', 'logg_L', 'logg_R', 'logg_Xc', 'logg_model']

REPO_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


# Unit conversions (must mirror src/frontend/ff_amp2.fypp)

def to_physical(d):
    return {'M': d[0] + 0.75,
            'Z': 10**(0.8*d[1] - 2.1),
            'Y': 0.1*d[2] + 0.22,
            'alpha': 2.0*d[3] + 1.0}


def to_scaled(M, Y, Z, alpha):
    return [M - 0.75, (math.log10(Z) + 2.1)/0.8, (Y - 0.22)/0.1, (alpha - 1.0)/2.0]


def key_of(d):
    return tuple(round(x, 9) for x in d)


def tc_str(T_c):
    """Canonical string for a list of candidate pre-MS T_c ('' = EM default)"""
    return ','.join(f'{t:.6g}' for t in T_c) if T_c else ''


class Prober:

    def __init__(self, args):

        self.exe = os.path.abspath(args.exe)
        self.stage = args.stage if args.stage != 'logg' else f'logg{args.log_g:g}'
        self.log_g = args.log_g if args.stage == 'logg' else None
        self.f_ov = args.f_ov
        self.max_age = args.max_age
        self.timeout = args.timeout
        self.threads = args.threads
        self.out = args.out
        self.log_dir = args.log_dir
        self.verbose = args.verbose
        self.T_c = tc_str(args.T_c)
        self.lock = threading.Lock()
        self.cache = {}

        if self.log_dir:
            os.makedirs(self.log_dir, exist_ok=True)

        if os.path.exists(self.out):
            with open(self.out, newline='') as f:
                for row in csv.DictReader(f):
                    if float(row['f_ov']) != self.f_ov or float(row['max_age']) != self.max_age:
                        continue
                    d = [float(row['d_'+a]) for a in AXES]
                    self.cache[(row['stage'], row.get('T_c') or '', key_of(d))] = row
            with open(self.out, newline='') as f:
                header = next(csv.reader(f), [])
            if header != FIELDS:
                sys.exit(f'{self.out} was written by an older version of this script; use a new --out file')
        else:
            with open(self.out, 'w', newline='') as f:
                csv.DictWriter(f, FIELDS).writeheader()

    def lookup(self, d, T_c=None):

        k = key_of(d)
        tc = self.T_c if T_c is None else T_c

        # A log g stage is only answered by itself; a zams-stage result
        # also answers a create-stage query

        if self.stage.startswith('logg'):
            return self.cache.get((self.stage, tc, k))

        row = self.cache.get(('zams', tc, k))
        if row is None and self.stage == 'create':
            row = self.cache.get(('create', tc, k))
        if row is not None and self.stage == 'create' and row['outcome'] == ZAMS_FAIL:
            row = dict(row, outcome=OK)

        return row

    def run(self, d, T_c=None):
        """Run (or fetch from cache) the probe at scaled point d, with candidate pre-MS
        central temperatures T_c (a tc_str; None = the prober's default); return the outcome"""

        tc = self.T_c if T_c is None else T_c

        row = self.lookup(d, tc)
        if row is not None:
            return row['outcome']

        p = to_physical(d)

        cmd = [self.exe, repr(p['M']), repr(p['Y']), repr(p['Z']), repr(p['alpha']),
               repr(self.f_ov), repr(self.max_age)]
        if tc:
            cmd.append(f'T_c={tc}')
        if self.log_g is not None:
            cmd.append(f'logg={self.log_g!r}')

        env = dict(os.environ, OMP_NUM_THREADS=str(self.threads))

        work_dir = tempfile.mkdtemp(prefix='em_probe_')

        t_start = time.time()

        stdout_lines = []
        t_code = ''
        create_ok = False
        T_c_used = tc.split(',')[0] if tc else 'default[1]'
        zams_props = {}
        logg_t_code = ''
        timed_out = []

        with open(os.path.join(work_dir, 'stderr'), 'w+') as err:

            proc = subprocess.Popen(cmd, cwd=work_dir, env=env, stdout=subprocess.PIPE,
                                    stderr=err, text=True, bufsize=1)

            def on_timeout():
                timed_out.append(True)
                proc.kill()

            timer = threading.Timer(self.timeout, on_timeout) if self.timeout > 0 else None
            if timer:
                timer.start()

            for line in proc.stdout:
                stdout_lines.append(line)
                if line.startswith('EM_PROBE CREATE_OK'):
                    create_ok = True
                    if self.stage == 'create':
                        proc.kill()
                        break
                elif line.startswith('EM_PROBE ZAMS_PROPS'):
                    zams_props = dict(zip(ZAMS_FIELDS, line.split()[2:]))
                elif line.startswith('EM_PROBE ZAMS_DONE'):
                    t_code = line.split()[2]
                elif line.startswith('EM_PROBE LOGG_PROPS'):
                    zams_props.update(zip(LOGG_FIELDS, line.split()[2:]))
                elif line.startswith('EM_PROBE LOGG_DONE'):
                    logg_t_code = line.split()[2]
                elif 'pre-MS model converged with T_c =' in line:
                    T_c_used = f'{float(line.split("=")[-1]):.6g}'

            proc.stdout.close()
            exit_code = proc.wait()

            if timer:
                timer.cancel()

            err.seek(0)
            stderr_text = err.read()

        runtime = time.time() - t_start

        # Classify

        if self.stage == 'create' and create_ok:
            outcome = OK
        elif timed_out:
            outcome = TIMEOUT
        elif not create_ok and exit_code < 0:
            outcome = CRASH
        elif not create_ok:
            outcome = CREATE_FAIL
        elif t_code == '':
            outcome = CRASH
        elif t_code != '0':
            outcome = ZAMS_FAIL
        elif self.log_g is not None and logg_t_code == '':
            outcome = CRASH
        elif self.log_g is not None and logg_t_code != '0':
            outcome = LOGG_FAIL
        else:
            outcome = OK

        message = ''
        if outcome != OK:
            m = re.findall(r'assertion .*|failed .*|ERROR.*|Program received signal.*', stderr_text)
            message = m[-1].strip() if m else (stderr_text.strip().splitlines() or [''])[-1].strip()
            if exit_code < 0 and outcome != TIMEOUT:
                message = f'killed by signal {-exit_code}; ' + message

        row = {'stage': self.stage, 'M': p['M'], 'Y': p['Y'], 'Z': p['Z'], 'alpha': p['alpha'],
               'f_ov': self.f_ov, 'max_age': self.max_age, 'T_c': tc, 'T_c_used': T_c_used if create_ok else '',
               'outcome': outcome, 't_code': t_code,
               'exit_code': exit_code, 'runtime': f'{runtime:.1f}', 'message': message[:200]}
        for a, x in zip(AXES, d):
            row['d_'+a] = repr(x)
        row.update(zams_props)
        row['logg_t_code'] = logg_t_code

        # Save logs

        if self.log_dir:
            tag = '_'.join(f'{x:.6f}' for x in d) + (f'_Tc{tc.replace(",", "-")}' if tc else '')
            with open(os.path.join(self.log_dir, f'{self.stage}_{tag}.log'), 'w') as f:
                f.write(' '.join(cmd) + '\n\n--- stdout ---\n' + ''.join(stdout_lines))
                f.write('\n--- stderr ---\n' + stderr_text)

        shutil.rmtree(work_dir, ignore_errors=True)

        with self.lock:
            self.cache[(self.stage, tc, key_of(d))] = row
            with open(self.out, 'a', newline='') as f:
                csv.DictWriter(f, FIELDS).writerow(row)
            if self.verbose:
                tc_note = f'  T_c={tc}' if tc else ''
                print(f'  [{outcome:>11s}] {fmt_point(d)}{tc_note}  ({runtime:.0f} s) {message}', flush=True)

        return outcome

    def run_many(self, points, jobs, T_c_list=None):
        """Run all points, for each entry of T_c_list (tc_strs; default the prober's T_c)"""
        points = list({key_of(d): d for d in points}.values())
        tasks = [(d, tc) for tc in (T_c_list or [self.T_c]) for d in points]
        with cf.ThreadPoolExecutor(max_workers=jobs) as ex:
            return list(ex.map(lambda t: self.run(*t), tasks))


def fmt_point(d):
    p = to_physical(d)
    return ('d=(' + ', '.join(f'{x:.5f}' for x in d) + ')  ' +
            f'M={p["M"]:.5f} Y={p["Y"]:.5f} Z={p["Z"]:.6f} alpha={p["alpha"]:.5f}')


def axis_offsets(width, n):
    if n < 2:
        return [0.0]
    return [-width + 2*width*i/(n-1) for i in range(n)]


def shifted(base, i, delta):
    d = list(base)
    d[i] = base[i] + delta
    return d


ENFORCE_BOUNDS = True


def in_bounds(d):
    return not ENFORCE_BOUNDS or all(0.0 <= x <= 1.0 for x in d)


def failure_extent(offsets, outcomes):
    """Return the contiguous failing interval of offsets containing 0"""
    i0 = min(range(len(offsets)), key=lambda i: abs(offsets[i]))
    if outcomes[i0] in (OK, None):
        return None
    lo = hi = i0
    while lo > 0 and outcomes[lo-1] not in (OK, None):
        lo -= 1
    while hi < len(offsets)-1 and outcomes[hi+1] not in (OK, None):
        hi += 1
    return (offsets[lo], offsets[hi], lo == 0, hi == len(offsets)-1)


# Modes

def do_point(prober, base, args):
    outcome = prober.run(base)
    print(f'\n{fmt_point(base)}\n  outcome: {outcome}')
    row = prober.lookup(base)
    if row and row['message']:
        print(f'  message: {row["message"]}')


def do_scan(prober, base, args):

    offsets = axis_offsets(args.width, args.n)

    points = []
    for i in range(4):
        for off in offsets:
            d = shifted(base, i, off)
            if in_bounds(d):
                points.append(d)

    print(f'Scanning {len({key_of(d) for d in points})} points (width ±{args.width} scaled units)...')
    prober.run_many(points, args.jobs)

    print('\nKey: ' + '  '.join(f'{s}={o}' for o, s in SYMBOL.items()) + '  (blank=out of bounds)')
    print(f'Offsets: {offsets[0]:+.4g} ... {offsets[-1]:+.4g} (step {offsets[1]-offsets[0]:.4g})\n')

    for i, a in enumerate(AXES):
        outcomes = []
        for off in offsets:
            d = shifted(base, i, off)
            outcomes.append(prober.run(d) if in_bounds(d) else None)
        line = ''.join(SYMBOL[o] if o else ' ' for o in outcomes)
        run = [o for o in outcomes if o is not None]
        n_fail = sum(o != OK for o in run)
        desc = f'{n_fail}/{len(run)} fail'

        # Report a failure interval only when the failures form a
        # single contiguous block through the base point; otherwise
        # the failures are scattered and an interval would mislead

        ext = failure_extent(offsets, outcomes)
        if ext is not None:
            lo, hi, lo_open, hi_open = ext
            n_block = sum(lo - 1e-12 <= off <= hi + 1e-12 for off in offsets)
            if n_block == n_fail:
                lo_phys = to_physical(shifted(base, i, lo))[a]
                hi_phys = to_physical(shifted(base, i, hi))[a]
                desc += (f'; contiguous: d in [{lo:+.4g}{"(or beyond)" if lo_open else ""}, '
                         f'{hi:+.4g}{"(or beyond)" if hi_open else ""}]  '
                         f'-> {a} in [{lo_phys:.6g}, {hi_phys:.6g}]')
            else:
                desc += '; scattered'
        elif n_fail > 0:
            desc += '; base point passes, scattered'
        print(f'  {a:>5s} |{line}|  {desc}')

    distinct = {key_of(d): d for d in points}.values()
    n_all = sum(prober.lookup(d)['outcome'] != OK for d in distinct)
    n_pts = len(distinct)
    print(f'\n  Overall: {n_all}/{n_pts} distinct points fail')


def do_edges(prober, base, args):

    base_outcome = prober.run(base)
    base_fail = base_outcome != OK
    print(f'Base point: {base_outcome}')

    def search(i, sign):

        # Coarse outward march until the outcome flips, then bisect

        step = args.width/args.n
        inner = 0.0
        outer = None
        k = 1
        while k*step <= args.width + 1e-12:
            d = shifted(base, i, sign*k*step)
            if not in_bounds(d):
                break
            if (prober.run(d) != OK) != base_fail:
                outer = k*step
                break
            inner = k*step
            k += 1

        if outer is None:
            return (i, sign, inner, None)

        while outer - inner > args.tol:
            mid = 0.5*(inner + outer)
            if (prober.run(shifted(base, i, sign*mid)) != OK) != base_fail:
                outer = mid
            else:
                inner = mid

        return (i, sign, inner, outer)

    tasks = [(i, s) for i in range(4) for s in (-1, +1)]

    with cf.ThreadPoolExecutor(max_workers=max(1, min(args.jobs, len(tasks)))) as ex:
        results = list(ex.map(lambda t: search(*t), tasks))

    what = 'failure' if base_fail else 'success'
    print(f'\nExtent of the {what} region around the base point (scaled offsets; '
          f'boundary bracketed to ±{args.tol:g}):\n')
    for i, sign, inner, outer in results:
        a = AXES[i]
        dirn = '+' if sign > 0 else '-'
        if outer is None:
            reach = shifted(base, i, sign*inner)
            print(f'  {a:>5s} {dirn}: no flip within {inner:.4g} '
                  f'(up to {a}={to_physical(reach)[a]:.6g})')
        else:
            b = shifted(base, i, sign*0.5*(inner+outer))
            print(f'  {a:>5s} {dirn}: flips between {inner:.4g} and {outer:.4g} '
                  f'(boundary at {a}≈{to_physical(b)[a]:.6g})')


def do_grid(prober, base, args):

    ix, iy = (AXES.index(a) for a in args.axes)
    offsets = axis_offsets(args.width, args.n)

    def pt(ox, oy):
        d = list(base)
        d[ix] += ox
        d[iy] += oy
        return d

    points = [pt(ox, oy) for oy in offsets for ox in offsets if in_bounds(pt(ox, oy))]
    print(f'Running {len(points)}-point grid in ({args.axes[0]}, {args.axes[1]})...')
    prober.run_many(points, args.jobs)

    print('\nKey: ' + '  '.join(f'{s}={o}' for o, s in SYMBOL.items()) + '  (blank=out of bounds)')
    print(f'x: d_{args.axes[0]} offset, y: d_{args.axes[1]} offset (base point marked by [ ])\n')
    for oy in reversed(offsets):
        cells = []
        for ox in offsets:
            d = pt(ox, oy)
            s = SYMBOL[prober.run(d)] if in_bounds(d) else ' '
            cells.append(f'[{s}]' if abs(ox) < 1e-12 and abs(oy) < 1e-12 else f' {s} ')
        print(f'  {oy:+8.4f} ' + ''.join(cells))
    ax, ay = args.axes
    print(f'\n  x: d_{ax} offset {offsets[0]:+.4g} ... {offsets[-1]:+.4g} -> '
          f'{ax} = {to_physical(pt(offsets[0], 0))[ax]:.6g} ... {to_physical(pt(offsets[-1], 0))[ax]:.6g}')
    print(f'  y: d_{ay} offset {offsets[0]:+.4g} ... {offsets[-1]:+.4g} -> '
          f'{ay} = {to_physical(pt(0, offsets[0]))[ay]:.6g} ... {to_physical(pt(0, offsets[-1]))[ay]:.6g}')


def do_tc(prober, base, args):

    # Points: a sweep along one axis (absolute scaled range), or the
    # 1-D scan neighbourhood of the base point

    rows = []
    if args.axis:
        i = AXES.index(args.axis)
        lo, hi = args.range
        for k in range(args.n):
            x = lo + (hi - lo)*k/(args.n - 1) if args.n > 1 else lo
            d = list(base)
            d[i] = x
            rows.append((f'{args.axis}={to_physical(d)[args.axis]:.5g}', d))
    else:
        for i, a in enumerate(AXES):
            for off in axis_offsets(args.width, args.n):
                d = shifted(base, i, off)
                if in_bounds(d):
                    rows.append((f'{a} {off:+.4f}', d))
        seen = set()
        rows = [r for r in rows if not (key_of(r[1]) in seen or seen.add(key_of(r[1])))]

    tcs = [tc_str([t]) if t > 0 else '' for t in args.tc_values]

    print(f'Running {len(rows)} points x {len(tcs)} T_c values...')
    prober.run_many([d for _, d in rows], args.jobs, tcs)

    print('\nKey: ' + '  '.join(f'{s}={o}' for o, s in SYMBOL.items()))
    print('Columns: pre-MS central temperature T_c (each tried on its own)\n')

    w = max(len(lbl) for lbl, _ in rows) + 2
    print(' '*w + ''.join(f'{t or "default":>8s}' for t in tcs) + '   any')
    n_ok = [0]*len(tcs)
    n_any = 0
    for lbl, d in rows:
        outs = [prober.run(d, tc) for tc in tcs]
        for j, o in enumerate(outs):
            n_ok[j] += (o == OK)
        any_ok = any(o == OK for o in outs)
        n_any += any_ok
        print(f'{lbl:<{w}s}' + ''.join(f'{SYMBOL[o]:>8s}' for o in outs) + f'{"yes" if any_ok else "NO":>6s}')
    print(' '*w + ''.join(f'{f"{k}/{len(rows)}":>8s}' for k in n_ok) + f'{f"{n_any}/{len(rows)}":>6s}')


def main():

    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('mode', choices=['point', 'scan', 'edges', 'grid', 'tc'])
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument('--scaled', nargs=4, type=float, metavar=('dM', 'dZ', 'dY', 'dalpha'),
                   help='base point in AMP scaled units (as passed to testffamp2)')
    g.add_argument('--physical', nargs=4, type=float, metavar=('M', 'Y', 'Z', 'alpha'),
                   help='base point in physical units (note order: M Y Z alpha)')
    ap.add_argument('--stage', choices=['create', 'zams', 'logg'], default='zams',
                    help='stop once the initial pre-MS model is built (~40%% cheaper), after reaching the ZAMS '
                         '(default), or after then evolving to log g = --log-g')
    ap.add_argument('--log-g', type=float, default=4.0, help='target log g for --stage logg (default 4.0)')
    ap.add_argument('--width', type=float, default=0.02,
                    help='half-width of the explored region, in scaled units (default 0.02)')
    ap.add_argument('-n', type=int, default=9,
                    help='points per axis for scan/grid; coarse steps per direction for edges (default 9)')
    ap.add_argument('--tol', type=float, default=1e-3, help='bisection tolerance for edges (scaled units)')
    ap.add_argument('--axes', nargs=2, choices=AXES, default=['M', 'Y'], help='axes for grid mode')
    ap.add_argument('--T_c', nargs='+', type=float, default=None,
                    help='candidate pre-MS central temperatures, tried in turn (default: EM\'s built-in list)')
    ap.add_argument('--tc-values', nargs='+', type=float,
                    default=[1e5, 1.5e5, 2e5, 3e5, 4e5, 5e5, 6e5, 8e5],
                    help='tc mode: pre-MS central temperatures to try, each on its own (0 = EM default list)')
    ap.add_argument('--axis', choices=AXES, default=None,
                    help='tc mode: sweep this axis over --range instead of using the scan neighbourhood')
    ap.add_argument('--range', nargs=2, type=float, default=[0.0, 1.0], metavar=('LO', 'HI'),
                    help='tc mode: scaled range for --axis (default 0 1)')
    ap.add_argument('-j', '--jobs', type=int, default=4, help='concurrent probe processes (default 4)')
    ap.add_argument('--threads', type=int, default=2, help='OMP_NUM_THREADS per probe (default 2)')
    ap.add_argument('--timeout', type=float, default=1800, help='per-run timeout in s (0 = none)')
    ap.add_argument('--f-ov', type=float, default=0.0, help='overshoot f_ov (ff_amp2 uses 0)')
    ap.add_argument('--max-age', type=float, default=1.5e10, help='max age (ff_amp2 uses 1.5e10)')
    ap.add_argument('--out', default='probe_results.csv', help='CSV results cache (default probe_results.csv)')
    ap.add_argument('--log-dir', default=None, help='if set, keep full stdout/stderr of each run here')
    ap.add_argument('--exe', default=os.path.join(REPO_DIR, 'bin', 'em_probe'), help='path to em_probe')
    ap.add_argument('-q', '--quiet', dest='verbose', action='store_false', help='no per-run progress lines')
    args = ap.parse_args()

    if not os.access(args.exe, os.X_OK):
        sys.exit(f'{args.exe} not found; build EM with FRONTENDS=yes')
    if 'MESA_DIR' not in os.environ:
        sys.exit('MESA_DIR must be set (MESA reads its data files at run time)')

    if args.scaled:
        base = list(args.scaled)
    else:
        base = to_scaled(*args.physical)

    print(f'Base point: {fmt_point(base)}')
    print(f'Stage: {args.stage}   f_ov={args.f_ov}   max_age={args.max_age:g}   cache: {args.out}')
    global ENFORCE_BOUNDS
    if not in_bounds(base):
        print('Warning: base point lies outside the AMP scaled range [0,1]; not clipping to it')
        ENFORCE_BOUNDS = False

    prober = Prober(args)

    {'point': do_point, 'scan': do_scan, 'edges': do_edges, 'grid': do_grid,
     'tc': do_tc}[args.mode](prober, base, args)


if __name__ == '__main__':
    main()

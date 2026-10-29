# -*- coding: utf-8 -*-
"""Run every probe against an economicspace checkout and compare its headline figures.

The notes are a measurement record: each figure belongs to the catalog release
it was taken on, and a re-run against a later release moves it (CLAUDE.md).
Until now the drift was found by hand and written into README. This runs each
probe, reads its headline figures back out of its output with the patterns in
`figures.json`, and prints them beside the recorded values. It also writes the
run to a dated JSON file outside this checkout, so two runs can be compared
later, and it says which release the checkout had installed.

It changes nothing in either repo. Each probe still reads only.

A probe whose self-check prints MISMATCH, or whose figure cannot be found in
its output, fails the run, and the exit status says so. The slow Lambert probe
(about ten minutes) runs only with --slow, and the probes that call public
APIs can be left out with --offline.

Usage, from this repo with economicspace beside it:
    py run_all.py
    py run_all.py --repo ../economicspace --slow
    py run_all.py --only spin covariance
    py run_all.py --offline
    py run_all.py --replay ~/.cache/starred-repos/snapshots/2026-10-29

Network responses are recorded through snapshots.py on a live run, and
--replay reads a recorded folder instead, so the API-backed figures can be
reproduced after the services have moved on.
"""

import argparse
import datetime
import json
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import snapshots  # noqa: E402
DEFAULT_OUT = os.path.join(os.path.expanduser("~"), ".cache", "starred-repos", "runs")


def load_spec(path=os.path.join(HERE, "figures.json")):
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def release_tag(repo):
    """The catalog release the checkout has installed, from its manifest."""
    path = os.path.join(repo, "asteroid_pipeline", "catalog_manifest.json")
    try:
        with open(path, encoding="utf-8") as fh:
            return json.load(fh).get("release_tag", "unknown")
    except (OSError, ValueError):
        return "unknown"


def run_probe(repo, script, args, timeout_s):
    """Run one probe with the checkout as its working directory; (output, seconds, status)."""
    cmd = [sys.executable, os.path.join(HERE, script)] + list(args)
    env = dict(os.environ, PYTHONIOENCODING="utf-8")
    start = time.time()
    try:
        done = subprocess.run(cmd, cwd=repo, capture_output=True, text=True,
                              encoding="utf-8", errors="replace", timeout=timeout_s,
                              stdin=subprocess.DEVNULL, env=env)
        out = done.stdout + done.stderr
        status = "ok" if done.returncode == 0 else f"exit {done.returncode}"
    except subprocess.TimeoutExpired as err:
        out = (err.stdout or "") if isinstance(err.stdout, str) else ""
        status = "timeout"
    if "MISMATCH" in out:
        status = "self-check failed"
    return out, time.time() - start, status


def extract(pattern, text):
    m = re.search(pattern, text)
    if not m:
        return None
    return float(m.group(1).replace(",", ""))


def main():
    ap = argparse.ArgumentParser(description="run every probe and compare its headline figures")
    ap.add_argument("--repo", default=os.path.join(os.path.dirname(HERE), "economicspace"),
                    help="economicspace checkout (default: beside this repo)")
    ap.add_argument("--only", nargs="*", help="probe keys from figures.json")
    ap.add_argument("--slow", action="store_true", help="include the slow probes")
    ap.add_argument("--offline", action="store_true", help="skip probes that call public APIs")
    ap.add_argument("--out", default=DEFAULT_OUT, help="where the dated JSON goes")
    ap.add_argument("--timeout", type=int, default=1800, help="seconds per probe")
    ap.add_argument("--keep-output", action="store_true",
                    help="also save each probe's full output beside the JSON")
    ap.add_argument("--replay", default=None,
                    help="a snapshots folder to replay network responses from (snapshots.py)")
    args = ap.parse_args()
    if args.replay:
        os.environ[snapshots.REPLAY_ENV] = os.path.abspath(args.replay)

    spec = load_spec()
    repo = os.path.abspath(args.repo)
    if not os.path.isfile(os.path.join(repo, "asteroid_pipeline", "asteroid_catalog.csv")):
        raise SystemExit(f"{repo} has no installed catalog; point --repo at an economicspace checkout")

    chosen = []
    for key, p in spec["probes"].items():
        if args.only and key not in args.only:
            continue
        if p.get("slow") and not args.slow and not args.only:
            continue
        if p.get("network") and args.offline:
            continue
        chosen.append(key)

    tag = release_tag(repo)
    stamp = datetime.datetime.now().strftime("%Y-%m-%dT%H%M%S")
    print(f"economicspace at {repo}, catalog release {tag}")
    print(snapshots.describe())
    print(f"probes: {', '.join(chosen)}\n")

    outputs, statuses, seconds = {}, {}, {}
    for key in chosen:
        p = spec["probes"][key]
        print(f"  running {p['script']} {' '.join(p['args'])} ...", flush=True)
        outputs[key], seconds[key], statuses[key] = run_probe(repo, p["script"], p["args"],
                                                              args.timeout)
        print(f"    {statuses[key]} in {seconds[key]:.0f} s")

    rows, failed = [], [k for k, s in statuses.items() if s != "ok"]
    for f in spec["figures"]:
        if f["probe"] not in outputs:
            continue
        now = extract(f["pattern"], outputs[f["probe"]])
        if now is None:
            failed.append(f"{f['probe']}: figure not found ({f['label']})")
        rows.append(dict(f, now=now))

    print(f"\n{'finding':7s} {'figure':58s} {'recorded':>10s} {'now':>10s} {'moved':>10s}  recorded on")
    for r in rows:
        moved = "" if r["now"] is None else f"{r['now'] - r['recorded']:+.4g}"
        now = "missing" if r["now"] is None else f"{r['now']:.6g}"
        print(f"{r['finding']:7s} {r['label'][:58]:58s} {r['recorded']:>10.6g} {now:>10s}"
              f" {moved:>10s}  {r['recorded_on']}")

    os.makedirs(args.out, exist_ok=True)
    path = os.path.join(args.out, f"run-{stamp}.json")
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        folder = snapshots.replay_folder() or snapshots.today_folder()
        try:
            manifest = snapshots.manifest_sha256(folder)
        except OSError:
            manifest = None
        json.dump({"run_at": stamp, "release_tag": tag, "repo": repo,
                   "network": {"mode": "replay" if snapshots.replay_folder() else "live",
                               "folder": folder, "manifest_sha256": manifest},
                   "probes": {k: {"status": statuses[k], "seconds": round(seconds[k], 1)}
                              for k in chosen},
                   "figures": [{k: r[k] for k in ("finding", "probe", "label", "recorded",
                                                  "recorded_on", "now")} for r in rows]},
                  fh, indent=1)
        fh.write("\n")
    if args.keep_output:
        for key, text in outputs.items():
            with open(os.path.join(args.out, f"run-{stamp}-{key}.txt"), "w",
                      encoding="utf-8", newline="\n") as fh:
                fh.write(text)
    print(f"\nwrote {path}")

    if failed:
        print("\nFAILED:")
        for f in failed:
            print(f"  {f}")
        sys.exit(1)


if __name__ == "__main__":
    main()

# -*- coding: utf-8 -*-
"""Record every network response a probe reads, and replay a recording later.

NHATS, Asterank, NEODyS and Horizons answer differently from one week to the
next, so a figure taken from them cannot be reproduced once the service has
moved on. The Toliou table solved this with a recorded sha256 (F12). This
does the same for every API call: each response is saved under a dated
folder with its URL, size and sha256 in a manifest, and a later run can replay
that folder instead of going to the network.

    live (default)    fetch, save to ~/.cache/starred-repos/snapshots/<date>/,
                      and add the response to that folder's manifest.json
    replay            STARRED_REPOS_REPLAY=<folder>: read every response from
                      the folder, check its sha256 against the manifest, and
                      fail rather than fetch if a URL was never recorded

The folder lives outside this checkout, because this repo sits in OneDrive and
some responses run to megabytes. A finding that rests on a recording names its
folder and the manifest's own sha256, so the recording can be shared and
checked.

Only the standard library is used.
"""

import datetime
import hashlib
import io
import json
import os
import re
import urllib.error
import urllib.request

ROOT = os.path.join(os.path.expanduser("~"), ".cache", "starred-repos", "snapshots")
REPLAY_ENV = "STARRED_REPOS_REPLAY"
USER_AGENT = "starred-repos probes (read-only; github.com/loggger101/starred-repos)"


class NotRecorded(LookupError):
    """Replay mode was asked for a URL its folder never recorded."""


class Tampered(ValueError):
    """A recorded response no longer has the sha256 its manifest names."""


def _sha256(data):
    return hashlib.sha256(data).hexdigest()


def _safe_name(url):
    name = re.sub(r"[^A-Za-z0-9._-]+", "_", url.split("://", 1)[-1])
    return name[:150] + "_" + _sha256(url.encode())[:12]


def _load_manifest(folder):
    path = os.path.join(folder, "manifest.json")
    if not os.path.exists(path):
        return {}
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def _save_manifest(folder, manifest):
    path = os.path.join(folder, "manifest.json")
    tmp = path + ".part"
    with open(tmp, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(manifest, fh, indent=1, sort_keys=True)
        fh.write("\n")
    os.replace(tmp, path)


def replay_folder():
    """The folder being replayed, or None for a live run."""
    return os.environ.get(REPLAY_ENV) or None


def today_folder(root=ROOT):
    return os.path.join(root, datetime.date.today().isoformat())


def fetch(url, timeout=180, root=ROOT, opener=None):
    """The response body for `url`, recorded on a live run and replayed otherwise.

    `opener` exists for the tests: a callable taking (url, timeout) and
    returning bytes, used in place of the network.
    """
    folder = replay_folder()
    if folder:
        manifest = _load_manifest(folder)
        entry = manifest.get(url)
        if entry is None:
            raise NotRecorded(f"{url} is not in the recording at {folder}")
        if entry.get("status"):
            # A recorded HTTP error (a 404 for a body the service does not
            # know) replays as the same error, so a "missing" count reproduces.
            raise urllib.error.HTTPError(url, entry["status"], "recorded", None, io.BytesIO())
        with open(os.path.join(folder, entry["file"]), "rb") as fh:
            data = fh.read()
        if _sha256(data) != entry["sha256"]:
            raise Tampered(f"{entry['file']} in {folder} does not match its recorded sha256")
        return data

    folder = today_folder(root)
    os.makedirs(folder, exist_ok=True)
    try:
        if opener is None:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=timeout) as fh:
                data = fh.read()
        else:
            data = opener(url, timeout)
    except urllib.error.HTTPError as err:
        manifest = _load_manifest(folder)
        manifest[url] = {"status": err.code, "file": None,
                         "fetched_at": datetime.datetime.now().isoformat(timespec="seconds")}
        _save_manifest(folder, manifest)
        raise

    name = _safe_name(url)
    with open(os.path.join(folder, name), "wb") as fh:
        fh.write(data)
    manifest = _load_manifest(folder)
    manifest[url] = {"file": name, "bytes": len(data), "sha256": _sha256(data),
                     "fetched_at": datetime.datetime.now().isoformat(timespec="seconds")}
    _save_manifest(folder, manifest)
    return data


def describe(root=ROOT):
    """One line naming where this run's responses come from or go to."""
    folder = replay_folder()
    if folder:
        return f"replaying network responses from {folder}"
    return f"recording network responses to {today_folder(root)}"


def manifest_sha256(folder):
    """The sha256 of a recording's manifest, which pins every response in it."""
    with open(os.path.join(folder, "manifest.json"), "rb") as fh:
        return _sha256(fh.read())

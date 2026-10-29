"""snapshots.py records responses with their sha256 and replays them, and refuses a tampered one."""

import io
import os
import tempfile
import unittest
import urllib.error
from unittest import mock

import tests._util  # noqa: F401  (puts the repo root on sys.path)
import snapshots


def fake_opener(responses):
    def opener(url, timeout):
        body = responses[url]
        if isinstance(body, int):
            raise urllib.error.HTTPError(url, body, "fake", None, io.BytesIO())
        return body
    return opener


class Snapshots(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = self.tmp.name
        self.env = mock.patch.dict(os.environ, {snapshots.REPLAY_ENV: ""})
        self.env.start()
        self.responses = {"https://example.test/a.json": b'{"x": 1}',
                          "https://example.test/missing": 404}
        opener = fake_opener(self.responses)
        snapshots.fetch("https://example.test/a.json", root=self.root, opener=opener)
        with self.assertRaises(urllib.error.HTTPError):
            snapshots.fetch("https://example.test/missing", root=self.root, opener=opener)
        self.folder = snapshots.today_folder(self.root)

    def tearDown(self):
        self.env.stop()
        self.tmp.cleanup()

    def replaying(self):
        return mock.patch.dict(os.environ, {snapshots.REPLAY_ENV: self.folder})

    def test_replay_returns_the_recorded_body(self):
        with self.replaying():
            self.assertEqual(snapshots.fetch("https://example.test/a.json"), b'{"x": 1}')

    def test_replay_reproduces_a_recorded_http_error(self):
        with self.replaying(), self.assertRaises(urllib.error.HTTPError) as err:
            snapshots.fetch("https://example.test/missing")
        self.assertEqual(err.exception.code, 404)

    def test_replay_refuses_an_unrecorded_url(self):
        with self.replaying(), self.assertRaises(snapshots.NotRecorded):
            snapshots.fetch("https://example.test/never")

    def test_replay_refuses_a_tampered_response(self):
        manifest = snapshots._load_manifest(self.folder)
        name = manifest["https://example.test/a.json"]["file"]
        with open(os.path.join(self.folder, name), "wb") as fh:
            fh.write(b'{"x": 2}')
        with self.replaying(), self.assertRaises(snapshots.Tampered):
            snapshots.fetch("https://example.test/a.json")

    def test_manifest_records_size_and_hash(self):
        entry = snapshots._load_manifest(self.folder)["https://example.test/a.json"]
        self.assertEqual(entry["bytes"], 8)
        self.assertEqual(entry["sha256"], snapshots._sha256(b'{"x": 1}'))
        self.assertEqual(len(snapshots.manifest_sha256(self.folder)), 64)


if __name__ == "__main__":
    unittest.main()

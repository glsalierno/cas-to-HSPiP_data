import unittest
from unittest.mock import patch

import requests

import pubchem_http as http


def response(status, retry_after=None):
    result = requests.Response()
    result.status_code = status
    result.url = "https://pubchem.ncbi.nlm.nih.gov/example"
    if retry_after is not None:
        result.headers["Retry-After"] = str(retry_after)
    return result


class HttpTests(unittest.TestCase):
    def setUp(self):
        self.patch = patch.object(http, "_blocked_until", 0)
        self.patch.start()
        self.addCleanup(self.patch.stop)

    def test_retry_429_and_503_with_timeout_and_backoff(self):
        with patch.object(http.time, "sleep") as sleep, patch.object(http.time, "monotonic", return_value=100), \
             patch.object(http.requests, "get", side_effect=[response(429, 3), response(503), response(200)]) as get:
            result = http.pubchem_get("https://example.test", timeout=30, context="test")
        self.assertEqual(result.status_code, 200)
        self.assertEqual(get.call_count, 3)
        self.assertTrue(all(call.kwargs == {"timeout": 30} for call in get.call_args_list))
        self.assertIn(((3,), {}), sleep.call_args_list)
        self.assertIn(((2,), {}), sleep.call_args_list)

    def test_retries_exhausted(self):
        with patch.object(http.time, "sleep"), patch.object(http.requests, "get", return_value=response(503)) as get:
            with self.assertRaises(requests.HTTPError):
                http.pubchem_get("https://example.test", timeout=30, context="test")
        self.assertEqual(get.call_count, 3)

    def test_403_does_not_retry(self):
        with patch.object(http.time, "sleep"), patch.object(http.requests, "get", return_value=response(403)) as get:
            with self.assertRaises(requests.HTTPError):
                http.pubchem_get("https://example.test", timeout=30, context="test")
        self.assertEqual(get.call_count, 1)

    def test_long_retry_after_stops_later_requests_too(self):
        with patch.object(http.time, "sleep"), patch.object(http.requests, "get", return_value=response(429, 90)) as get:
            with self.assertRaises(requests.HTTPError):
                http.pubchem_get("https://example.test/one", timeout=30, context="test")
            with self.assertRaisesRegex(requests.RequestException, "cooldown"):
                http.pubchem_get("https://example.test/two", timeout=30, context="test")
        self.assertEqual(get.call_count, 1)

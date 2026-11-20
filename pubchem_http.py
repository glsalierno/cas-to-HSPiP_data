"""Small synchronous PubChem request policy, shared by this repo's CLI."""
from datetime import datetime, timezone
from email.utils import parsedate_to_datetime
import logging
import time

import requests

LOG = logging.getLogger(__name__)
TRANSIENT_STATUS = {429, 503}
MAX_ATTEMPTS = 3
MAX_RETRY_AFTER_SECONDS = 30
_blocked_until = 0.0


def retry_after_seconds(header, default):
    if header:
        try:
            return max(0, int(header))
        except ValueError:
            try:
                deadline = parsedate_to_datetime(header)
            except (ValueError, TypeError, OverflowError):
                return default
            if deadline.tzinfo is not None:
                return max(0, (deadline - datetime.now(timezone.utc)).total_seconds())
    return default


def pubchem_get(url, *, timeout, context):
    """Return a successful Response; raise on a final HTTP or network error."""
    global _blocked_until
    for attempt in range(1, MAX_ATTEMPTS + 1):
        remaining = _blocked_until - time.monotonic()
        if remaining > MAX_RETRY_AFTER_SECONDS:
            raise requests.RequestException("PubChem cooldown active; no request was sent")
        if remaining > 0:
            time.sleep(remaining)
        # Allow at most four request starts per second in this sequential CLI.
        time.sleep(0.25)
        response = requests.get(url, timeout=timeout)
        if response.status_code not in TRANSIENT_STATUS:
            response.raise_for_status()
            return response
        delay = max(2 ** (attempt - 1), retry_after_seconds(response.headers.get("Retry-After"), 0))
        _blocked_until = time.monotonic() + delay
        if attempt == MAX_ATTEMPTS:
            LOG.warning("%s: PubChem HTTP %d after %d attempts", context, response.status_code, attempt)
            response.raise_for_status()
        if delay > MAX_RETRY_AFTER_SECONDS:
            LOG.warning("%s: PubChem asks for a %.1f-second wait; stopping", context, delay)
            response.raise_for_status()
        LOG.warning("%s: PubChem HTTP %d (attempt %d/%d); retrying in %.1fs",
                    context, response.status_code, attempt, MAX_ATTEMPTS, delay)

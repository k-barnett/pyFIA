"""
Southwest Region MOG helpers (AZ / NM; R ``region %in% "southwest"``).

- :mod:`fia_mog.southwest.core` — ``HABTYPCD1`` → ecological response unit (ERU).
- :mod:`fia_mog.southwest.evaluate` — Table 9 / 19 MOG vector (used by the engine).

``southwest_eru`` remains importable from :mod:`fia_mog.crosswalk` for a single
crosswalk surface. The package ``__init__`` only pulls in :mod:`core` so that
``crosswalk`` can import ``southwest_eru`` without a circular import through
``engine``.
"""

from __future__ import annotations

from .core import southwest_eru

__all__ = [
    "southwest_eru",
]

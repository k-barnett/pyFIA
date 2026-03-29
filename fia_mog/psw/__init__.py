"""
Pacific Southwest (California / R ``region %in% "pacific southwest"``) MOG logic.

Table 12 old growth and Table 19 maturity are implemented in ``psw.evaluate``.
NWFP overlap for **white fir** uses the same boundary as the Pacific Northwest
(``utility_NWFPboundary``), populated on :attr:`fia_mog.engine.ConditionContext.pnw_inside_nwfp`.
"""

from __future__ import annotations

__all__: list[str] = []

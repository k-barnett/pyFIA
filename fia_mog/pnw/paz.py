"""
Plant association zone (PAZ) raster value → group label.

Mirrors ``FUNCTION_mapMOG.R`` lines 2987–3003.
"""

from __future__ import annotations

import numpy as np


def paz_value_to_group(paz_val: int | float | None) -> str | None:
    if paz_val is None or (isinstance(paz_val, float) and np.isnan(paz_val)):
        return None
    v = int(paz_val)
    if v in {56, 57, 58, 59}:
        return "white fir - grand fir"
    if v in {40, 41, 42, 43}:
        return "douglas fir"
    if v in {25, 27, 29}:
        return "lodgepole pine"
    if v in {87, 88, 89}:
        return "silver fir"
    if v in {30, 31, 33, 34}:
        return "ponderosa pine"
    if v in {61, 62, 63, 64}:
        return "subalpine fir"
    if v in {82, 83, 84}:
        return "western hemlock"
    if v in {91, 92, 93, 94}:
        return "mountain hemlock"
    if v in {21, 22}:
        return "juniper"
    if v == 18:
        return "oak woodland"
    if v == 79:
        return "port orford cedar"
    if v in {48, 49}:
        return "redwood"
    if v in {71, 72, 73}:
        return "shasta red fir"
    if v in {46, 47}:
        return "sitka spruce"
    if v in {51, 52}:
        return "tanoak"
    if v in {37, 38, 39}:
        return "jeffrey pine - knobcone pine"
    return None

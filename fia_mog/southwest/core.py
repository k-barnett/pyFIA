"""
Southwest Region (R ``region %in% "southwest"``) crosswalk helpers — AZ / NM, ADFOR 3xx.

ERU mapping from condition ``HABTYPCD1`` mirrors ``FUNCTION_mapMOG.R`` (~2187–2207).
"""

from __future__ import annotations

import numpy as np


def southwest_eru(habtypcd1: int | float | None) -> str | None:
    """
    Map HABTYPCD1 to Ecological Response Unit (ERU) for the Southwest rules.

    Ports the chained ``if`` assignments in ``FUNCTION_mapMOG.R`` (lines ~2187–2207).
    The script assigns ``pinyon juniper sagebrush`` and then re-assigns the **same**
    code vector to ``pinyon juniper grass``; the effective ERU for those codes is
    **grass** (Table 9 SDI threshold), not sagebrush (which would use QMD).
    """

    if habtypcd1 is None or np.isnan(habtypcd1):
        return None
    h = int(habtypcd1)

    if h in {
        415,
        435,
        604,
        1100,
        3060,
        3080,
        3090,
        3110,
        3111,
        3112,
        3200,
        3201,
        3202,
        3203,
        3231,
        3240,
        3300,
        3301,
        3310,
        3320,
        3350,
        3370,
        3999,
        4060,
        4061,
        4062,
        4151,
        4152,
        4300,
        4310,
        4320,
        4330,
        4340,
        4350,
        4351,
        4360,
        4999,
        26005,
        240300,
    }:
        return "spruce-fir forest"
    if h in {
        1010,
        1011,
        1012,
        1020,
        1030,
        1070,
        1080,
        1081,
        1110,
        1111,
        1120,
        1150,
        1160,
        1231,
        1999,
        6010,
        6060,
        6070,
        6071,
        6080,
        6130,
        12320,
        12333,
    }:
        return "mixed conifer with aspen"
    if h in {238040, 238310}:
        return "bristlecone pine"
    if h in {
        1021,
        1022,
        1040,
        1041,
        1042,
        1050,
        1051,
        1052,
        1053,
        1054,
        1060,
        1090,
        1140,
        1141,
        1203,
        1213,
        1239,
        1241,
        6090,
        11130,
        12140,
        12141,
        12142,
        12143,
        12330,
        12331,
        12332,
        12340,
        12341,
        12350,
        12360,
        12361,
        12362,
        12380,
        12420,
        12430,
        12999,
        238300,
    }:
        return "mixed conifer - frequent fire"
    if h in {
        11030,
        11031,
        11032,
        11033,
        11035,
        11090,
        11091,
        11092,
        11093,
        11210,
        11211,
        11212,
        11213,
        11214,
        11215,
        11216,
        11320,
        11330,
        11340,
        11341,
        11350,
        11380,
        11390,
        11391,
        11392,
        11400,
        11460,
        11500,
        11999,
    }:
        return "ponderosa pine forest"
    if h in {
        11034,
        11220,
        11360,
        11361,
        11370,
        11410,
        11411,
        11420,
        11430,
        11440,
        32010,
        32030,
        32999,
        33010,
        33020,
        33030,
    }:
        return "ponderosa pine - evergreen oak"
    if h in {
        3102,
        204400,
        230030,
        230040,
        230041,
        230042,
        230999,
        231010,
        232070,
        233010,
        233030,
        233040,
        233041,
        233042,
        233050,
    }:
        return "pinyon juniper evergreen shrub"
    if h in {202500, 204320, 204330, 204500, 232020, 232330, 233330}:
        return "pinyon juniper (persistent)"
    if h in {20404, 204050, 204321, 2040303}:
        return "pinyon juniper deciduous shrub"
    # R: same c(...) as "pinyon juniper sagebrush" then immediately overwritten by "grass"
    if h in {
        20406,
        20410,
        20411,
        20431,
        23204,
        204021,
        204022,
        204023,
        204024,
        204300,
        204350,
        204370,
        204999,
        231020,
        232030,
        232999,
        233020,
        233021,
        233022,
        233999,
        9000042,
    }:
        return "pinyon juniper grass"
    if h in {
        20140,
        201010,
        201011,
        201020,
        201040,
        201331,
        201332,
        201333,
        201340,
        201350,
        201400,
        201410,
        201999,
        202320,
        202321,
        202330,
        202331,
        202999,
        231021,
        231030,
        231040,
        231050,
        231999,
        9000043,
    }:
        return "juniper grass"
    if h in {
        3101,
        204360,
        232050,
        232060,
        630010,
        630030,
        630040,
        630043,
        630050,
        2040301,
        2040302,
    }:
        return "madrean pinyon-oak"
    if h in {
        31999,
        610010,
        610020,
        620010,
        620020,
        620021,
        620030,
        620999,
        630020,
        630041,
        630042,
        632999,
        650010,
        650999,
    }:
        return "madrean encinal woodland"
    if h == 640999:
        return "gambel oak shrubland"
    if h in {201420, 201430, 210999}:
        return "semi-desert grassland"
    if h == 11470:
        return "ponderosa pine/willow"
    if h in {1130, 620040}:
        return "arizona walnut"
    if h == 104:
        return "rio grande cottonwood/shrub"
    if h == 103:
        return "narrowleaf cottonwood-spruce, narrowleaf cottonwood/shrub"
    if h == 3:
        return "upper montane conifer/willow"

    return "other"

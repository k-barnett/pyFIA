"""
Habitat veg-code → letter rules for Northern Region / western Montana zone (east of plot = west of divide in R naming).

Mirrors the R `northern.subregion == "western Montana zone"` branch.
"""

from __future__ import annotations

WEST_MT_RULES: list[tuple[str, frozenset[str]]] = [
    (
        "A",
        frozenset(
            """
            PIFL/AGSP PIFL/FIED PIFL/FESC PIFL/JUCO PIPO/AND PIPO/AGSP PIPO/FEID PIPI/FEID
            PIPO/FESC PIPO/PUTR PIPO/SYAL PIPO/PRVI PIPO/PRIV PIPO/SHCA PIPO/PHMA
            PSME/AGSP PSME/FEID PSME/FESC PSME/SYAL PSME/CARU PSME/ARUV
            """.split()
        ),
    ),
    (
        "B",
        frozenset(
            """
            PSME/VACA PSME/PHMA PSME/CARU PSME/SMST PSME/VAGL PSME/ARUV PSME/SYAL PSME/PIPO
            PSME/SPBE PSME/SYOR ABGR/SPBE ABGR/PHMA ABGR/COOC
            """.split()
        ),
    ),
    (
        "C",
        frozenset(
            """
            PSME/VAGL PSME/XETE PSME/LIBO PSME/CARU PSME/JUCO PSME/ARCO ABGR/XETE ABGR/COOC
            ABGR/VAGL ABLA/CARU
            """.split()
        ),
    ),
    (
        "D",
        frozenset(
            """
            ABGR/VAGL ABGR/ASCR ABGR/ASCA ABGR/MEFE ABGR/TABR ABGR/CLUN ABGR/ARNU ABGR/XETE
            ABGR/PHMA ABGR/SETR THPL/CLUN THPL/ARNU THPL/MEFE THPL/XETE THPL/TABR THPL/ASCA
            THPL/GYDR TSHE/GYDR TSHE/CLUN TSHE/ARNU TSHE/MEFE THSE/XETE TSHE/ASCA
            """.split()
        ),
    ),
    (
        "E",
        frozenset(
            """
            PIAB/CLUN PIBR/CLUN PIEN/CLUN PIGL/CLUN PIMA/CLUN PIPU/CLUN PIRU/CLUN PISI/CLUN PIOM/CLUN
            PIAB/VACA PIBR/VACA PIEN/VACA PIGL/VACA PIMA/VACA PIPU/VACA PIRU/VACA PISI/VACA PIOM/VACA
            PIAB/SEST PIBR/SEST PIEN/SEST PIGL/SEST PIMA/SEST PIPU/SEST PIRU/SEST PISI/SEST PIOM/SEST
            PIAB/PSME PIBR/PSME PIEN/PSME PIGL/PSME PIMA/PSME PIPU/PSME PIRU/PSME PISI/PSME PIOM/PSME
            ABLA/CLUN ABLA/ARNU ABLA/VACA ABLA/XETE ABLA/MEFE ABLA/LIBO ABLA/COOC ABLA/LUHI ABLA/VASC
            TSME/STAM TSME/MEFE TSME/LUHI TSME/XETE TSME/CLUN ABLA/ALSI ABLA/LUHI ABLA/MEFE
            """.split()
        ),
    ),
    (
        "F",
        frozenset(
            """
            PIAB/EQAR PIBR/EQAR PIEN/EQAR PIGL/EQAR PIMA/EQAR PIPU/EQAR PIRU/EQAR PISI/EQAR PIOM/EQAR
            PIAB/GATR PIBR/GATR PIEN/GATR PIGL/GATR PIMA/GATR PIPU/GATR PIRU/GATR PISI/GATR PIOM/GATR
            PIAB/SMST PIBR/SMST PIEN/SMST PIGL/SMST PIMA/SMST PIPU/SMST PIRU/SMST PISI/SMST PIOM/SMST
            THPL/ALFI THPL/ADPE THPL/OPHO ABLA/OPHO ABLA/GATR ABLA/CACA ABLA/STAM ABLA/MEFE ABLA/LICA
            ABLA/VACA ABLA/LEGL TSME/STAM TSME/LUHI TSME/MEFE
            """.split()
        ),
    ),
    ("G", frozenset("PSME/LIBO ABGR/LIBO ABGR/UBO PSME/VAGL".split())),
    (
        "H",
        frozenset(
            """
            PIAB/PHMA PIBR/PHMA PIEN/PHMA PIGL/PHMA PIMA/PHMA PIPU/PHMA PIRU/PHMA PISI/PHMA PIOM/PHMA
            PIAB/VACA PIBR/VACA PIEN/VACA PIGL/VACA PIMA/VACA PIPU/VACA PIRU/VACA PISI/VACA PIOM/VACA
            ABLA/VACA ABLA/LIBO ABLA/VASC ABLA/XETE ABLA/VAGL ABLA/COOC ABLA/LUHI TSME/XETE TSME/LUHI
            TSME/VAGL ABLA/THOC ABLA/ARCO ABLA/CAGE ABLA/GAGE ABLA/PSME ABLA/RIMO ABAL/RIMO
            PICO/PUTR PICO/VACA PICO/XETE PICO/LIBO PICO/VASC PICO/CARU
            """.split()
        ),
    ),
    ("I", frozenset("ABLA/VASC ABLA/PIAL ABLA/LUHI TSME/LUHI TSME/VASC TSME/MEFE".split())),
    ("J", frozenset("PIAL/ABLA PIAL LALY/ABLA".split())),
]


def west_mt_habitat_letters(veg_code: str) -> list[str]:
    return [letter for letter, codes in WEST_MT_RULES if veg_code in codes]

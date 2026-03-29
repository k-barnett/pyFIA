"""
Habitat veg-code → letter rules for Northern Region / eastern Montana zone.

Mirrors the R `northern.subregion == "eastern Montana zone"` branch.
"""

from __future__ import annotations

EAST_MT_RULES: list[tuple[str, frozenset[str]]] = [
    (
        "A",
        frozenset(
            """
            PIFL PIFL/AGSP PIFL/FEID PIFL/FESC PIFL/JUCO PIPO PIPO/JUSC PIPO/MARE PIPO/AND PIPO/AGSP
            PIPO/FEID PIPO/FESC PIPO/SYOC JUSC/PSSP JUSC/PSSPS PIPO/PUTR PIPO/SYAL PIPO/PRVI PIPO/SHCA
            PIPO/PHMA PSME/AGSP PSME/FEID PSME/FESC PSME/SYAL PSME/CARU PSME/ARUV PSME/JUCO PSME/ARCO
            PSME/SYOR PSME/JUSC
            """.split()
        ),
    ),
    (
        "B",
        frozenset(
            """
            PSME/PHMA PSME/SMST PSME/CARU PSME/ARUV PSME/PIPO PSME/CAGE PSME/SPBE PICO/PUTR
            """.split()
        ),
    ),
    (
        "C",
        frozenset(
            """
            PSME/MARE PSME/VACA PSME/PHMA PSME/VAGL PSME/SYAL
            PIAB/PHMA PIBR/PHMA PIEN/PHMA PIGL/PHMA PIMA/PHMA PIPU/PHMA PIRU/PHMA PISI/PHMA PIOM/PHMA
            """.split()
        ),
    ),
    (
        "D",
        frozenset(
            """
            PSME/VAGL PSME/XETE PSME/LIBO PSME/SYAL PSME/CARU
            PIAB/LIBO PIBR/LIBO PIEN/LIBO PIGL/LIBO PIMA/LIBO PIPU/LIBO PIRU/LIBO PISI/LIBO PIOM/LIBO
            PIAB/SMST PIBR/SMST PIEN/SMST PIGL/SMST PIMA/SMST PIPU/SMST PIRU/SMST PISI/SMST PIOM/SMST
            ABLA/LIBO ABLA/XETE ABLA/VASC PICO/LIBO
            """.split()
        ),
    ),
    (
        "E",
        frozenset(
            """
            PIAB/EQAR PIBR/EQAR PIEN/EQAR PIGL/EQAR PIMA/EQAR PIPU/EQAR PIRU/EQAR PISI/EQAR PIOM/EQAR
            PIAB/CLUN PIBR/CLUN PIEN/CLUN PIGL/CLUN PIMA/CLUN PIPU/CLUN PIRU/CLUN PISI/CLUN PIOM/CLUN
            PIAB/VACA PIBR/VACA PIEN/VACA PIGL/VACA PIMA/VACA PIPU/VACA PIRU/VACA PISI/VACA PIOM/VACA
            PIAB/GATR PIBR/GATR PIEN/GATR PIGL/GATR PIMA/GATR PIPU/GATR PIRU/GATR PISI/GATR PIOM/GATR
            ABLA/CLUN ABLA/ARNU ABLA/VACA ABLA/XETE ABLA/MEFE ABLA/GATR ALBL/GATR ALBL/VASC ABLA/CACA
            ABLA/LEGL
            """.split()
        ),
    ),
    (
        "F",
        frozenset(
            """
            PIAB/VACA PIBR/VACA PIEN/VACA PIGL/VACA PIMA/VACA PIPU/VACA PIRU/VACA PISI/VACA PIOM/VACA
            ABLA/SYAL ABLA/VACA ABLA/XETE ABLA/VAGL ABLA/VASC TSME/XETE ABLA/CARU ABLA/THOC PICO/VACA
            PICO/VASC PICO/CARU PICO/JUCO
            """.split()
        ),
    ),
    ("G", frozenset("ABLA/MEFE ABLA/ALSI".split())),
    (
        "H",
        frozenset(
            """
            PIAB/SEST PIBR/SEST PIEN/SEST PIGL/SEST PIMA/SEST PIPU/SEST PIRU/SEST PISI/SEST PIOM/SEST
            PICA/JUCO ABLA/JUCO ABLA/CARU ABLA/CLPS PIAB/PSME PIBR/PSME PIEN/PSME PIGL/PSME PIMA/PSME
            PIPU/PSME PIRU/PSME PISI/PSME PIOM/PSME ABLA/ARCO ABLA/CAGE ABLA/PSME
            """.split()
        ),
    ),
    (
        "I",
        frozenset(
            """
            TSME/MEFE ABLA/RIMO ABLA/VASC ABLA'LUHI ABLA/MEFE TSME/LUHI TSME/VASC TSME/XETE
            """.split()
        ),
    ),
    ("J", frozenset("PIAL-ABLA LALY-ABLA PIAL".split())),
]


def east_mt_habitat_letters(veg_code: str) -> list[str]:
    return [letter for letter, codes in EAST_MT_RULES if veg_code in codes]

"""
Habitat veg-code → letter rules for Northern Region / Idaho zone (Green et al. 1992 Appendix A).

Mirrors the R `northern.subregion == "northern Idaho zone"` branch.
"""

from __future__ import annotations

# (letter, frozenset of veg.code strings) — order matches R (multiple letters per code possible).
IDAHO_RULES: list[tuple[str, frozenset[str]]] = [
    (
        "A",
        frozenset(
            """
            PIPO PIPO/AGSP PIPO/FEID PIPO/FESC PIPO/PURT PIPO/SYAL PIPO/BERE
            PSME/AGSP PSME/FEID PSME/FESC
            """.split()
        ),
    ),
    (
        "B",
        frozenset(
            """
            PIPO/PRIV PIPO/SHCA PSME/PHMA PSME/SMST PSME/CARU PSME/SYAL PSME/VAGL
            VAGL/ARUV VAGL/XETE PSME/VACA PSME/AGSP PSME/CARU PSME/ARUV PSME/PIPO
            PSME/CAGE PSME/SPBE PSME/JUCO PSME/ARCO PSME/SYOR
            PIAB/PHMA PIBR/PHMA PIEN/PHMA PIGL/PHMA PIMA/PHMA PIPU/PHMA PIRU/PHMA PISI/PHMA
            """.split()
        ),
    ),
    (
        "C",
        frozenset(
            """
            ABGR/SETR ABGR/ASCA ABGR/MEFE ABGR/CLUN ABLA/CLUN ABLA/ARNU ABGR/PHMA
            """.split()
        ),
    ),
    ("C1", frozenset("ABGR/ASCA ABGR/TABR ABGR/CLUN".split())),
    (
        "D",
        frozenset(
            """
            PSME/LIBO PSME/SYAL PSME/CARU PSME/VAGL ABGR/LIBO ABGR/XETE ABGR/COOC
            ABGR/VAGL ABGR/CLUN
            """.split()
        ),
    ),
    ("E", frozenset("ABGR/PHMA ABGR/COOC ABGR/SPBE".split())),
    ("F", frozenset("THPL/OPHO THPL/ATFI THPL/ADPE".split())),
    (
        "G",
        frozenset(
            """
            THPL/CLUN THPL/ARNU THPL/GYDR THPL/ASCA THPL/MEFE THPL/XETE TSHE/GYDR
            TSHE/ASCA TSHE/ARNU TSHE/MEFE TSHE/CLUN TSHE/XETE
            """.split()
        ),
    ),
    ("G1", frozenset("THPL/CLUN THPL/TABR THPL/ASCA".split())),
    (
        "H",
        frozenset(
            """
            PIAB/EQAR PIBR/EQAR PIEN/EQAR PIGL/EQAR PIMA/EQAR PIPU/EQAR PIRU/EQAR PISI/EQAR
            PIAB/GATR PIBR/GATR PIEN/GATR PIGL/GATR PIMA/GATR PIPU/GATR PIRU/GATR PISI/GATR
            ABLA/OPHO ABLA/GATR ABLA/STAM ABLA/METE ABLA/LICA ABLA/CACA ABLA/VACA ABLA/LEGL
            TSME/STAM TSME/LUHI TSME/MEFE
            """.split()
        ),
    ),
    (
        "I",
        frozenset(
            """
            PIAB/CLUN PIBR/CLUN PIEN/CLUN PIGL/CLUN PIMA/CLUN PIPU/CLUN PIRU/CLUN PISI/CLUN
            PIAB/LIBO PIBR/LIBO PIEN/LIBO PIGL/LIBO PIMA/LIBO PIPU/LIBO PIRU/LIBO PISI/LIBO
            TSME/MEFE ABLA/CLUN ABLA/ARNU ABLA/VACA ABLA/XETE ABLA/MEFE ABLA/LIBO ABLA/COOC
            ABLA/LUHI ABLA/VASC TSME/CLUN TSME/XETE
            """.split()
        ),
    ),
    (
        "J",
        frozenset(
            """
            PIAB/CLUN PIBR/CLUN PIEN/CLUN PIGL/CLUN PIMA/CLUN PIPU/CLUN PIRU/CLUN PISI/CLUN
            ABLA/VACA ABLA/XETE ABLA/VAGL ABLA/VASC ABLA/COOC ABLA/LUHI TSME/XETE TSME/LUHI
            TSME/VASC ABLA/THOC ABLA/CARU ABLA/CLPS ABLA/ARCO ABLA/CAGE ABLA/PSME
            """.split()
        ),
    ),
    (
        "K",
        frozenset(
            """
            ABLA/RIMO ABLA-PIAL ABLA/VASC PIAL/VASC ABLA/LUHI TSME/LUHI ABLA/MEFE TSME/VASC
            TSME/MEFE PIAL-ABLA LALY-ABLA PIAL PICO/VACA PICO/XETE PICO/LIBO PICO/VASC PICO/CARU
            """.split()
        ),
    ),
]


def idaho_habitat_letters(veg_code: str) -> list[str]:
    return [letter for letter, codes in IDAHO_RULES if veg_code in codes]

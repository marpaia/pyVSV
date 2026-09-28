"""The reviews whose reference lists define the construction-paper bibliography."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class SeedReview:
    """One index review. ``pmcid`` is empty when the review is not in PMC."""

    pmid: str
    pmcid: str
    label: str


REVIEWS: tuple[SeedReview, ...] = (
    SeedReview("15121047", "", "Lichty 2004"),
    SeedReview("16977404", "PMC7079905", "Clarke 2006"),
    SeedReview("21987744", "PMC3218670", "Geisbert 2011"),
    SeedReview("27395563", "PMC5220644", "Monath 2016"),
    SeedReview("29757706", "PMC6183239", "Suder 2018"),
    SeedReview("31368826", "PMC6816421", "Fathi 2019"),
    SeedReview("34578125", "PMC8470541", "Liu 2021"),
)

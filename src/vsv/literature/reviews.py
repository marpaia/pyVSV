"""The reviews whose reference lists define the construction-paper bibliography."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class Review:
    """One index review. ``pmcid`` is empty when the review is not in PMC."""

    pmid: str
    pmcid: str
    label: str


REVIEWS: tuple[Review, ...] = (
    Review("15121047", "", "Lichty 2004"),
    Review("16977404", "PMC7079905", "Clarke 2006"),
    Review("21987744", "PMC3218670", "Geisbert 2011"),
    Review("27395563", "PMC5220644", "Monath 2016"),
    Review("29757706", "PMC6183239", "Suder 2018"),
    Review("31368826", "PMC6816421", "Fathi 2019"),
    Review("34578125", "PMC8470541", "Liu 2021"),
)

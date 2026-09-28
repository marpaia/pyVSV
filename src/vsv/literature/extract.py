"""Quote-backed architecture fields from a paper's abstract and, when present, its methods.

A field stays ``unknown`` when no sentence states it. Disagreeing sentences become a conflict
and the field stays ``unknown``. ΔG together with a named foreign antigen in the same sentence
is recorded as a replacement of the glycoprotein. A locus is recorded only from a stated position.
"""

from __future__ import annotations

import re
from collections.abc import Callable, Iterable

from vsv.design.evidence import Assay, Host, Outcome, Phenotype, Stage
from vsv.literature.classify import classify
from vsv.literature.fetch import clean_text
from vsv.literature.query import VECTOR_KINDS
from vsv.literature.records import UNKNOWN, Claim, Extraction, Paper, Quoted

_SENTENCE = re.compile(r"(?<=[.!?])\s+(?=[A-Z\"(])")
_DELTA_G = re.compile(
    r"vsv(?:\W{0,3})?(?:Δ|∆|delta)\W{0,3}g\b|\bg-deleted\b|glycoprotein-deleted",
    re.IGNORECASE,
)
_REPLACE = re.compile(
    r"replac\w* (?:the )?(?:vsv |viral )?(?:glycoprotein|\bg gene\b)|"
    r"in place of (?:the )?(?:vsv )?(?:glycoprotein|\bg gene\b)|"
    r"substitut\w* for (?:the )?(?:vsv )?glycoprotein",
    re.IGNORECASE,
)
_DELETE = re.compile(
    r"delet\w* (?:the )?(?:vsv |viral )?glycoprotein",
    re.IGNORECASE,
)
_MATRIX_POLYMERASE = re.compile(
    r"between the matrix(?: protein)? and (?:the )?(?:large )?(?:polymerase|\bl gene\b)",
    re.IGNORECASE,
)
_INSERT_BETWEEN = (
    (
        re.compile(
            r"between the nucleoprotein and (?:the )?phosphoprotein|"
            r"between (?:the )?n and p genes",
            re.IGNORECASE,
        ),
        "Nucleoprotein",
    ),
    (
        re.compile(
            r"between the phosphoprotein and (?:the )?matrix|"
            r"between (?:the )?p and m genes",
            re.IGNORECASE,
        ),
        "Phosphoprotein",
    ),
    (
        re.compile(
            r"between the matrix(?: protein)? and (?:the )?(?:glycoprotein|\bg gene\b)|"
            r"between (?:the )?m and g genes",
            re.IGNORECASE,
        ),
        "MatrixProtein",
    ),
    (
        re.compile(
            r"between the glycoprotein and (?:the )?(?:polymerase|\bl gene\b)|"
            r"between (?:the )?g and l genes",
            re.IGNORECASE,
        ),
        "Glycoprotein",
    ),
)
_AFTER = re.compile(
    r"(?:after|downstream of) the "
    r"(leader|nucleoprotein|phosphoprotein|matrix(?: protein)?|glycoprotein|"
    r"(?:large )?polymerase)",
    re.IGNORECASE,
)
_AFTER_GENE = {
    "leader": "leader",
    "nucleoprotein": "Nucleoprotein",
    "phosphoprotein": "Phosphoprotein",
    "matrix": "MatrixProtein",
    "matrix protein": "MatrixProtein",
    "glycoprotein": "Glycoprotein",
    "polymerase": "LargePolymerase",
    "large polymerase": "LargePolymerase",
}
_POSITION = re.compile(r"\bposition ([1-6])\b", re.IGNORECASE)
_POSITION_CONTEXT = re.compile(r"\b(?:gene|insert|transcription|proximal)\b", re.IGNORECASE)
_PROXIMAL = re.compile(r"3'[- ]proximal|first (?:gene )?position", re.IGNORECASE)
_INDIANA = re.compile(r"\bindiana\b", re.IGNORECASE)
_NEW_JERSEY = re.compile(r"\bnew jersey\b", re.IGNORECASE)
_SINGLE = re.compile(
    r"single[- ]cycle|replication[- ]incompetent|replication[- ]defective",
    re.IGNORECASE,
)
_COMPETENT = re.compile(r"replication[- ]competent", re.IGNORECASE)
_CLINICAL = re.compile(
    r"phase (?:1|i|2|ii|3|iii)\b|clinical trial|healthy adults",
    re.IGNORECASE,
)
_NEGATION = re.compile(r"\b(?:no|not|without|failed to|did not)\b", re.IGNORECASE)
_CHALLENGE = re.compile(
    r"\b(?:lethal challenge|challenged|protected from|survival)\b", re.IGNORECASE
)
_NEUTRAL = re.compile(r"\bneutraliz", re.IGNORECASE)
_IMMUNE = re.compile(r"\b(?:immunogen|antibod|seroconvert)", re.IGNORECASE)
_REACTO = re.compile(r"\breactogen", re.IGNORECASE)
_ADVERSE = re.compile(r"\b(?:adverse events?|hearing loss)\b", re.IGNORECASE)
_HOSTS: tuple[tuple[str, re.Pattern[str]], ...] = (
    (Host.NONHUMAN_PRIMATE.value, re.compile(r"non-?human primates?|\bmacaques?\b", re.IGNORECASE)),
    (Host.HAMSTER.value, re.compile(r"\bhamsters?\b", re.IGNORECASE)),
    (Host.MOUSE.value, re.compile(r"\bmice\b|\bmouse\b|\bmurine\b", re.IGNORECASE)),
    (
        Host.HUMAN.value,
        re.compile(r"\b(?:humans?|patients?|volunteers?|healthy adults)\b", re.IGNORECASE),
    ),
    (Host.CELL_LINE.value, re.compile(r"\bcell lines?\b|\bvero\b|\bbhk-?21\b", re.IGNORECASE)),
)
_ANIMAL = frozenset({Host.MOUSE.value, Host.HAMSTER.value, Host.NONHUMAN_PRIMATE.value})

# Specific names before broader ones. Reporter payloads are dropped when a pathogen slug is present.
_SLUGS: tuple[tuple[re.Pattern[str], str, str], ...] = (
    (
        re.compile(
            r"\bzebov[- ]?gpc?\b|\bebola virus glycoprotein\b|\bebolavirus glycoprotein\b",
            re.IGNORECASE,
        ),
        "filovirus-glycoprotein",
        "Ebola virus glycoprotein",
    ),
    (
        re.compile(r"\bsudv[- ]?gp\b|\bsudan virus glycoprotein\b", re.IGNORECASE),
        "sudan-ebolavirus-glycoprotein",
        "Sudan virus glycoprotein",
    ),
    (
        re.compile(r"\bbundibugyo\b.{0,40}\bglycoprotein\b|\bbdbv[- ]?gp\b", re.IGNORECASE),
        "bundibugyo-ebolavirus-glycoprotein",
        "Bundibugyo virus glycoprotein",
    ),
    (
        re.compile(r"\breston\b.{0,40}\bglycoprotein\b|\brestv[- ]?gp\b", re.IGNORECASE),
        "reston-ebolavirus-glycoprotein",
        "Reston virus glycoprotein",
    ),
    (
        re.compile(
            r"\bta[iï] forest\b.{0,40}\bglycoprotein\b|\btafv[- ]?gp\b",
            re.IGNORECASE,
        ),
        "taiforest-ebolavirus-glycoprotein",
        "Taï Forest virus glycoprotein",
    ),
    (
        re.compile(r"\bmarv[- ]?gp\b|\bmarburg\b.{0,40}\bglycoprotein\b", re.IGNORECASE),
        "marburgvirus-glycoprotein",
        "Marburg virus glycoprotein",
    ),
    (
        re.compile(r"\blasv[- ]?gpc\b|\blassa\b.{0,40}\b(?:gpc|glycoprotein)\b", re.IGNORECASE),
        "arenavirus-gpc",
        "Lassa virus glycoprotein complex",
    ),
    (
        re.compile(r"\bnipah\b.{0,40}\bfusion\b|\bfusion protein\b.{0,40}\bnipah\b", re.IGNORECASE),
        "henipavirus-fusion",
        "Nipah fusion protein",
    ),
    (
        re.compile(
            r"\bnipah\b.{0,40}\bglycoprotein\b|\bglycoprotein\b.{0,40}\bnipah\b",
            re.IGNORECASE,
        ),
        "henipavirus-glycoprotein",
        "Nipah glycoprotein",
    ),
    (
        re.compile(r"\bhendra\b.{0,40}\bfusion\b", re.IGNORECASE),
        "henipavirus-fusion",
        "Hendra fusion protein",
    ),
    (
        re.compile(r"\bhendra\b.{0,40}\bglycoprotein\b", re.IGNORECASE),
        "henipavirus-glycoprotein",
        "Hendra glycoprotein",
    ),
    (
        re.compile(
            r"\b(?:sars-cov-2|sars-cov2)\b.{0,50}\bspike\b|\bspike\b.{0,50}\b(?:sars-cov-2|sars-cov2)\b",
            re.IGNORECASE,
        ),
        "sarbecovirus-spike",
        "SARS-CoV-2 spike",
    ),
    (
        re.compile(r"\binfluenza\b.{0,40}\bhemagglutinin\b", re.IGNORECASE),
        "influenza-hemagglutinin",
        "influenza hemagglutinin",
    ),
    (
        re.compile(r"\bhiv-?1?\b.{0,40}\b(?:envelope|env)\b", re.IGNORECASE),
        "hiv-envelope",
        "HIV envelope",
    ),
    (
        re.compile(
            r"\b(?:rsv|respiratory syncytial virus)\b.{0,40}\bfusion\b",
            re.IGNORECASE,
        ),
        "pneumovirus-fusion",
        "RSV fusion protein",
    ),
    (
        re.compile(
            r"\bcrimean[- ]congo\b.{0,40}\bglycoprotein\b|\bcchfv\b.{0,30}\bglycoprotein\b",
            re.IGNORECASE,
        ),
        "nairovirus-glycoprotein",
        "Crimean-Congo hemorrhagic fever virus glycoprotein",
    ),
    (
        re.compile(
            r"\bzika\b.{0,40}\b(?:envelope|e protein)\b|\benvelope\b.{0,40}\bzika\b",
            re.IGNORECASE,
        ),
        "zika-envelope",
        "Zika envelope",
    ),
    (
        re.compile(
            r"\b(?:hantavirus|andes virus|hantaan)\b.{0,40}\bglycoprotein\b",
            re.IGNORECASE,
        ),
        "hantavirus-glycoprotein",
        "hantavirus glycoprotein",
    ),
    (
        re.compile(r"\bchikungunya\b.{0,40}\b(?:envelope|structural protein)\b", re.IGNORECASE),
        "chikungunya-envelope",
        "chikungunya envelope",
    ),
    (
        re.compile(r"\brabies\b.{0,40}\bglycoprotein\b", re.IGNORECASE),
        "rabies-glycoprotein",
        "rabies glycoprotein",
    ),
    (re.compile(r"\b(?:egfp|gfp)\b", re.IGNORECASE), "reporter-gfp", "GFP"),
    (re.compile(r"\bluciferase\b", re.IGNORECASE), "reporter-luciferase", "luciferase"),
)
_PATHOGENS: tuple[tuple[re.Pattern[str], str], ...] = (
    (re.compile(r"\bebola\b|\bzebov\b", re.IGNORECASE), "Ebola"),
    (re.compile(r"\bmarburg\b|\bmarv\b", re.IGNORECASE), "Marburg"),
    (re.compile(r"\blassa\b|\blasv\b", re.IGNORECASE), "Lassa"),
    (re.compile(r"\bnipah\b", re.IGNORECASE), "Nipah"),
    (re.compile(r"\bhendra\b", re.IGNORECASE), "Hendra"),
    (re.compile(r"\bsars-cov-2\b|\bsars-cov2\b", re.IGNORECASE), "SARS-CoV-2"),
    (re.compile(r"\binfluenza\b", re.IGNORECASE), "influenza"),
    (re.compile(r"\bhiv\b", re.IGNORECASE), "HIV"),
    (re.compile(r"\bzika\b", re.IGNORECASE), "Zika"),
    (re.compile(r"\bchikungunya\b", re.IGNORECASE), "chikungunya"),
    (re.compile(r"\brabies\b", re.IGNORECASE), "rabies"),
    (re.compile(r"\bhantavirus\b|\bandes virus\b", re.IGNORECASE), "hantavirus"),
)
_REPORTER_SLUGS = frozenset({"reporter-gfp", "reporter-luciferase"})


def sentences(text: str) -> list[str]:
    """Split ``text`` into sentences. Short fragments are kept."""
    cleaned = clean_text(text)
    if not cleaned:
        return []
    return [part.strip() for part in _SENTENCE.split(cleaned) if part.strip()]


def describe(paper: Paper, methods: str = "") -> Extraction:
    """Classify ``paper`` and extract architecture fields.

    ``methods`` is open-access methods text. It affects fields, not the class. The class is
    read from the title, abstract, publication type, and MeSH terms.
    """
    classification = classify(paper)
    source = "pmid-only"
    parts = [paper.title, paper.abstract]
    if paper.title or paper.abstract:
        source = "abstract"
    if methods.strip():
        parts.append(methods)
        source = "fulltext"
    body = sentences("\n".join(parts))
    if classification.kind not in VECTOR_KINDS:
        return _empty(
            paper.pmid, classification.kind, classification.review, classification.status, source
        )
    conflicts: list[str] = []
    backbone = _exclusive(body, _backbone, "backbone", conflicts)
    replication = _exclusive(body, _replication, "replication", conflicts)
    activity, target = _edit(body, conflicts)
    name, slug = _antigen(body, conflicts)
    locus = _exclusive(body, _locus, "locus", conflicts)
    claims = _claims(body)
    stage = _stage(body, claims)
    return Extraction(
        paper.pmid,
        classification.kind,
        classification.review,
        classification.status,
        backbone,
        replication,
        activity,
        target,
        name,
        slug,
        locus,
        stage,
        tuple(claims),
        tuple(conflicts),
        source,
    )


def _empty(pmid: str, kind: str, review: bool, status: str, source: str) -> Extraction:
    return Extraction(
        pmid,
        kind,
        review,
        status,
        UNKNOWN,
        UNKNOWN,
        UNKNOWN,
        UNKNOWN,
        UNKNOWN,
        UNKNOWN,
        UNKNOWN,
        UNKNOWN,
        (),
        (),
        source,
    )


def _exclusive(
    body: list[str],
    reader: Callable[[str], Quoted | None],
    field: str,
    conflicts: list[str],
) -> Quoted:
    found = [
        quoted for sentence in body if (quoted := reader(sentence)) is not None and quoted.known
    ]
    return _agree(found, field, conflicts)


def _agree(found: list[Quoted], field: str, conflicts: list[str]) -> Quoted:
    values = {item.value for item in found}
    if not values:
        return UNKNOWN
    if len(values) > 1:
        conflicts.append(field)
        return UNKNOWN
    return found[0]


def _backbone(sentence: str) -> Quoted | None:
    indiana = _INDIANA.search(sentence) is not None
    jersey = _NEW_JERSEY.search(sentence) is not None
    if indiana and jersey:
        return None
    if indiana:
        return Quoted("indiana", sentence)
    if jersey:
        return Quoted("new-jersey", sentence)
    return None


def _replication(sentence: str) -> Quoted | None:
    single = _SINGLE.search(sentence) is not None
    competent = _COMPETENT.search(sentence) is not None
    if single and competent:
        return None
    if single:
        return Quoted("single-cycle", sentence)
    if competent:
        return Quoted("replication-competent", sentence)
    return None


def _edit(body: list[str], conflicts: list[str]) -> tuple[Quoted, Quoted]:
    found: list[tuple[str, str, str]] = []
    for sentence in body:
        decision = _edit_sentence(sentence)
        if decision is not None:
            found.append(decision)
    if not found:
        for sentence in body:
            if _DELTA_G.search(sentence) and _payload(sentence):
                found.append(("replace-module", "Glycoprotein", sentence))
    activities = {item[0] for item in found}
    targets = {item[1] for item in found}
    if len(activities) > 1 or len(targets) > 1:
        conflicts.append("edit")
        return UNKNOWN, UNKNOWN
    if not found:
        return UNKNOWN, UNKNOWN
    activity, target, sentence = found[0]
    return Quoted(activity, sentence), Quoted(target, sentence)


def _edit_sentence(sentence: str) -> tuple[str, str, str] | None:
    deleted = _DELETE.search(sentence) is not None or _DELTA_G.search(sentence) is not None
    if _MATRIX_POLYMERASE.search(sentence) and deleted:
        return ("replace-module", "Glycoprotein", sentence)
    if _REPLACE.search(sentence):
        return ("replace-module", "Glycoprotein", sentence)
    inserted = _insert_target(sentence)
    if inserted is not None:
        return ("insert-module", inserted, sentence)
    if _DELETE.search(sentence):
        return ("delete-module", "Glycoprotein", sentence)
    return None


def _insert_target(sentence: str) -> str | None:
    for pattern, target in _INSERT_BETWEEN:
        if pattern.search(sentence):
            return target
    match = _AFTER.search(sentence)
    if match is None:
        return None
    return _AFTER_GENE.get(match.group(1).lower())


def _locus(sentence: str) -> Quoted | None:
    if _POSITION_CONTEXT.search(sentence) is None and _PROXIMAL.search(sentence) is None:
        return None
    if _PROXIMAL.search(sentence):
        return Quoted("0", sentence)
    match = _POSITION.search(sentence)
    if match is None:
        return None
    return Quoted(str(int(match.group(1)) - 1), sentence)


def _antigen(body: list[str], conflicts: list[str]) -> tuple[Quoted, Quoted]:
    slugs: list[tuple[str, str, str]] = []
    names: list[tuple[str, str]] = []
    for sentence in body:
        for pattern, slug, name in _SLUGS:
            if pattern.search(sentence):
                slugs.append((slug, name, sentence))
        for pattern, name in _PATHOGENS:
            if pattern.search(sentence):
                names.append((name, sentence))
    pathogen_slugs = [item for item in slugs if item[0] not in _REPORTER_SLUGS]
    chosen = pathogen_slugs or slugs
    unique = {item[0] for item in chosen}
    if len(unique) > 1:
        conflicts.append("antigen")
        return UNKNOWN, UNKNOWN
    if chosen:
        slug, name, sentence = chosen[0]
        return Quoted(name, sentence), Quoted(slug, sentence)
    unique_names = {item[0] for item in names}
    if len(unique_names) > 1:
        conflicts.append("antigen")
        return UNKNOWN, UNKNOWN
    if names:
        name, sentence = names[0]
        return Quoted(name, sentence), UNKNOWN
    return UNKNOWN, UNKNOWN


def _payload(sentence: str) -> bool:
    if any(pattern.search(sentence) for pattern, _slug, _name in _SLUGS):
        return True
    return any(pattern.search(sentence) for pattern, _name in _PATHOGENS)


def _claims(body: Iterable[str]) -> list[Claim]:
    claims: list[Claim] = []
    seen: set[tuple[str, str, str]] = set()
    for sentence in body:
        phenotype = _phenotype(sentence)
        if phenotype is None:
            continue
        kind, assay, outcome = phenotype
        for host, pattern in _HOSTS:
            if pattern.search(sentence) is None:
                continue
            key = (kind, host, assay)
            if key in seen:
                continue
            seen.add(key)
            claims.append(Claim(kind, host, assay, outcome, sentence))
    return claims


def _phenotype(sentence: str) -> tuple[str, str, str] | None:
    negated = _NEGATION.search(sentence) is not None
    outcome = Outcome.NOT_REPORTED.value if negated else Outcome.REPORTED.value
    if _ADVERSE.search(sentence):
        assay = Assay.REACTOGENICITY.value if _REACTO.search(sentence) else Assay.SAFETY.value
        harm = Outcome.ADVERSE.value if not negated else outcome
        return Phenotype.SAFETY.value, assay, harm
    if _REACTO.search(sentence):
        return Phenotype.SAFETY.value, Assay.REACTOGENICITY.value, outcome
    if _CHALLENGE.search(sentence):
        return Phenotype.PROTECTION.value, Assay.LETHAL_CHALLENGE.value, outcome
    if _NEUTRAL.search(sentence):
        return Phenotype.IMMUNOGENICITY.value, Assay.NEUTRALIZATION.value, outcome
    if _IMMUNE.search(sentence):
        return Phenotype.IMMUNOGENICITY.value, Assay.IMMUNOGENICITY.value, outcome
    return None


def _stage(body: list[str], claims: list[Claim]) -> Quoted:
    for sentence in body:
        if _CLINICAL.search(sentence):
            return Quoted(Stage.CLINICAL.value, sentence)
    for claim in claims:
        if claim.host in _ANIMAL and claim.outcome == Outcome.REPORTED.value:
            return Quoted(Stage.PRECLINICAL.value, claim.quote)
    return UNKNOWN

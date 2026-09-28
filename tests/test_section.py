"""Article text split into quotable blocks."""

from __future__ import annotations

from vsv.literature.documents.jats import section_xml
from vsv.literature.documents.text import section_text

_ARTICLE = """
<article>
  <front><article-meta>
    <abstract><title>Abstract</title><p>Short abstract.</p></abstract>
  </article-meta></front>
  <body>
    <sec><title>Results</title>
      <sec><title>Henipaviruses</title><p>Nipah was inserted.</p></sec>
      <fig><caption><p>Figure 1. The genome.</p></caption></fig>
    </sec>
    <sec sec-type="ack"><title>Acknowledgements</title><p>Thanks.</p></sec>
    <sec><title>Methods</title>
      <table-wrap>
        <caption><p>Table 1. Constructs.</p></caption>
        <table>
          <thead><tr><th>Construct</th><th>Backbone</th></tr></thead>
          <tbody>
            <tr><td>rVSV-SUDV</td><td>Indiana</td></tr>
            <tr><td>rVSV-BDBV</td><td>Indiana</td></tr>
          </tbody>
        </table>
      </table-wrap>
      <p>The insert was AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA and then stopped.</p>
    </sec>
  </body>
  <back><ref-list><ref><p>A cited paper.</p></ref></ref-list></back>
</article>
"""


def test_xml_sections_keep_tables_and_drop_the_reference_list() -> None:
    blocks = section_xml(_ARTICLE)
    assert [(block.section, block.kind, block.text) for block in blocks] == [
        ("Abstract", "paragraph", "Short abstract."),
        ("Results / Henipaviruses", "paragraph", "Nipah was inserted."),
        ("Results", "caption", "Figure 1. The genome."),
        ("Methods", "caption", "Table 1. Constructs."),
        ("Methods", "table-row", "Construct: rVSV-SUDV | Backbone: Indiana"),
        ("Methods", "table-row", "Construct: rVSV-BDBV | Backbone: Indiana"),
        (
            "Methods",
            "paragraph",
            "The insert was [sequence omitted] and then stopped.",
        ),
    ]
    assert [block.id for block in blocks] == [f"b{index:04d}" for index in range(1, 8)]
    assert "Thanks." not in " ".join(block.text for block in blocks)
    assert "A cited paper." not in " ".join(block.text for block in blocks)


def test_plain_text_uses_headings_and_column_lines() -> None:
    text = "\n".join(
        [
            "Results",
            "",
            "Mice were protected.",
            "",
            "Construct    Backbone",
            "rVSV-SUDV    Indiana",
        ]
    )
    blocks = section_text(text)
    assert [(block.section, block.kind, block.text) for block in blocks] == [
        ("Results", "paragraph", "Mice were protected."),
        ("Results", "table-row", "Construct | Backbone"),
        ("Results", "table-row", "rVSV-SUDV | Indiana"),
    ]

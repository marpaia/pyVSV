# Construction-paper census

`index.json` is the bibliography. It is the union of PubMed identifiers cited by seven reviews: Lichty 2004 (PMID 15121047), Clarke 2006 (PMC7079905), Geisbert 2011 (PMC3218670), Monath 2016 (PMC5220644), Suder 2018 (PMC6183239), Fathi 2019 (PMC6816421), and Liu 2021 (PMC8470541). Reference lists come from Europe PMC. A PMCID on a row comes from the PMC id converter, not from the frozen corpus flag.

`manifest.jsonl` is one download attempt per cited PMID. `status` is `xml-body`, `front-matter-only`, `no-pmc`, or `error`. An `xml-body` file is `data/corpus/fulltext/{PMCID}.xml`. A body is usable when its text is at least 2,000 characters. Front matter is not saved.

```sh
VSV_NCBI_EMAIL=name@example.edu vsv census index
VSV_NCBI_EMAIL=name@example.edu vsv census fetch
vsv census section
```

`vsv census section` writes `sections/{pmid}.json` for each `xml-body` row. A block is one paragraph, caption, or table row. Those files are gitignored because they contain the article.

`VSV_NCBI_EMAIL` is the contact address NCBI asks for on E-utilities requests. The commands sleep between requests. Re-running `vsv census fetch` leaves an existing XML file in place when its body still passes the length test.

`acquired/`, `sections/`, and `drafts/` are local caches and are gitignored. They are not part of this bibliography yet.

## Precedent scope

`precedents` compares three arrangements: a glycoprotein replacement, one extra gene immediately after VSV G, and one extra gene immediately after the leader. A plan with an unknown transcribed locus, or with any number of non-viral genes other than one, is not compared.

An empty result means one of the following.

- The candidate is not one of those three arrangements.
- A transcribed locus is unknown. `rvsv_wnv_e` and `rvsv_hcv_ce1e2` keep VSV G, and the saved construction paper does not state where the extra gene sits. On `rvsv_ch_wnv_e` the West Nile envelope locus is unknown.
- The genome has more than one non-viral module. `rvsv_niv_n` and `rvsv_n_zika_prme` add a gene to the filovirus glycoprotein replacement. `rvsv_ch_hiv_gag`, `rvsv_ch_hsv2_gd`, and `rvsv_nj_lcrv` each encode a replacement glycoprotein and one extra gene.
- No other curated plan in the families below has that arrangement. When a phenotype or host filter is set, an empty result means no recorded observation matches. Follow-up papers that use a genome already in the library are not all entered as observations.

The construction papers recorded here cover the Indiana chassis, reporter inserts, glycoprotein replacements, and extra transcription units for the antigens and effectors named on the plans. Those families include filovirus glycoproteins (Zaire, Sudan, Bundibugyo, and Bombali), Marburg, Lassa, Nipah F, G, and N, Hendra, SARS-CoV, SARS-CoV-2, MERS, chikungunya, SFTSV, hantavirus, CCHF, rabies, HIV envelope and gag, influenza, RSV, Zika, dengue-2, HSV-2 gD, CRPV L1, HCV glycoproteins, HCV core/E1/E2, West Nile envelope, Yersinia LcrV, GM-CSF, IFN-β, IL-4, and thymidine kinase. Chandipura G or New Jersey G encoded in place of Indiana G is its own plan when the paper states that swap. CAT, leader GFP, and leader GM-CSF are in the library. The G-kept chikungunya genome is `rvsv_gl_chikv_e`.

Some saved results sit on the plan and are not phenotype protection. A protection filter does not return them.

- Rabbit papilloma prevention by `rvsv_gl_crpv_l1` is unmodeled. It is not a lethal challenge.
- The IFN-β tumor result is unmodeled. The regression sentence in the IL-4 and thymidine-kinase paper does not name either virus, so it is not an observation.
- Survival after a prime with one vector and a boost with another is not protection by either genome alone. The saved regimens are Chandipura HSV-2 gD primed and Indiana gD boosted, LcrV primed from the leader or from between G and L and boosted with the New Jersey vector, and Indiana West Nile envelope primed and Chandipura envelope boosted.

A blend of separately built vectors is not one genome. Bundibugyo Group 4, and a mix of Sudan, Zaire, and Marburg vectors, are blends. A pseudotype whose entry protein is supplied outside the genome is not one of these plans. A matrix residue, a polymerase methyltransferase substitution, or a gene-order shuffle is not stored, because a viral-gene module keeps identity and locus only. An XhoI–NheI plasmid site is not an insert after G unless the same paper states the junction. The G-kept HCV E1/E2 clones do not state that junction.

`vsv precedents` prints the plan, phenotype, host, assay, and citation. The note is not in that table. A result that depends on a second vector is not entered as protection.

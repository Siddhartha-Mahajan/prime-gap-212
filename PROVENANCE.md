# Authorship, provenance, and verification status

This file records how the repository artifact was produced. It is intended to prevent the use of vague phrases such as "AI-assisted" from understating the model's role.

## Attribution

| Role | Attribution |
| --- | --- |
| Generator and repository-artifact author | OpenAI GPT-5.6 Sol (`gpt-5.6-sol`), used in the ChatGPT web interface with Pro reasoning effort |
| Later review and computational verification | Codex |
| Human prompter, curator, and releaser | Siddhartha Mahajan |
| Research resources | Lossfunk, headed by Paras Chopra |

Official model documentation: <https://developers.openai.com/api/docs/models/gpt-5.6-sol>

## What Siddhartha Mahajan did

Siddhartha Mahajan supplied the bounded-prime-gaps problem, pointed GPT-5.6 Sol in the ChatGPT web interface to the relevant papers and project resources, and gave brief high-level directions, including asking it to pursue the strongest reduction it could support. He selected the resulting package for public release and later asked Codex to inspect and test it.

Mahajan did **not** independently derive, line-check, or verify the mathematical proof. He did not write the manuscript or verification code. His questions to Codex and his execution of this release workflow are not independent mathematical verification.

## What GPT-5.6 Sol did

GPT-5.6 Sol, running in the ChatGPT web interface with Pro reasoning effort, conducted the literature-guided investigation, developed the claimed support and partition argument, searched for the finite sieve function, generated the mathematical exposition, produced the certificate and verification programs, and assembled the manuscript package. Codex was used only afterward to review the package and rerun computational checks.

The model is named as the author of this GitHub artifact at Mahajan's request because that best describes who generated it.

## Computational checks

The release contains two implementations of the finite calculation:

- a directed-rounding IEEE binary64 verifier; and
- a fixed-scale, checked 256-bit integer-interval verifier.

During a later review session, Codex reran the integer verifier on Apple ARM64 using Clang and pinned Boost 1.92.0. It returned

```text
1.0001949701565048417353464247383560
    <= 45 J_box(F) / I(F) <=
1.0001949701568054616041293775800977
```

and reported that the published enclosure and strict threshold were verified. These calculations support reproducibility of the finite certificate. They do not constitute independent human review of the analytic reduction from Stadlmann's theorems to the certificate.

## Current status

No independent human mathematician has yet verified the analytic proof. Until such review occurs, this repository should be described as an **AI-generated candidate mathematical result with reproducible computational certificates**, not as an independently verified or peer-reviewed theorem.

The repository includes this human-readable account but does not currently include a complete native export of the Codex conversation.

Canonical repository: <https://github.com/Siddhartha-Mahajan/prime-gap-212>

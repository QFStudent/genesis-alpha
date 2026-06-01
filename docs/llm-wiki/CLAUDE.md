# LLM Wiki — Quant Research Schema

This is the schema document for the genesis-alpha quant research wiki. You are the wiki maintainer. Obsidian is the IDE; you are the programmer; the wiki is the codebase.

**Your job:** Read, cross-reference, summarize, and file. The human curates sources and asks questions. You do the bookkeeping that makes the knowledge base useful over time.

---

## Purpose

Accumulate, synthesize, and cross-reference knowledge from academic papers, articles, and notes on topics that feed into genesis-alpha strategy development:

- Alpha signals and factors (momentum, carry, value, short-term reversal, quality, etc.)
- Market microstructure (market impact, liquidity, adverse selection, execution costs)
- Portfolio construction and risk management (optimization, risk parity, factor exposure)
- Backtesting methodology (bias types, multiple testing, deflated Sharpe, walk-forward)
- Futures-specific topics (roll, continuous construction, basis, term structure)

The goal is a compounding knowledge base — not re-derived answers every time, but a structured, maintained wiki that gets richer with every source ingested and every question filed.

---

## Directory Layout

```
llm-wiki/
├── CLAUDE.md              # This file — the schema and workflow instructions
├── index.md               # Content catalog (updated on every ingest)
├── log.md                 # Append-only activity log
├── raw/                   # Source documents — immutable, never modified by LLM
│   └── (PDFs, markdown clips, transcripts, notes)
└── wiki/                  # LLM-maintained pages — you own this layer
    ├── papers/            # One page per paper or article
    ├── concepts/          # Alpha signals, factors, risk metrics, microstructure
    ├── strategies/        # Strategy overviews and multi-paper syntheses
    └── entities/          # Instruments, venues, data vendors, datasets
```

External sources live at `../references/` — treat them the same as `raw/`.

---

## Page Frontmatter Convention

Every wiki page starts with YAML frontmatter:

```yaml
---
type: paper | concept | strategy | entity | synthesis
tags: [tag1, tag2]          # lowercase, hyphenated; use the tag list below
sources: [AuthorYYYY-slug]  # slugs of papers this page draws on
updated: YYYY-MM-DD
---
```

Use `[[PageName]]` for internal Obsidian wikilinks. For cross-directory links: `[[papers/Author2023]]`, `[[concepts/MarketImpact]]`.

### Approved Tags

**Asset class:** `equities`, `futures`, `fx`, `fixed-income`, `crypto`, `multi-asset`  
**Signal type:** `momentum`, `short-term-reversal`, `carry`, `value`, `quality`, `low-vol`, `event-driven`  
**Topic:** `market-impact`, `liquidity`, `execution`, `portfolio-construction`, `risk`, `factor-model`, `backtest-methodology`, `data`  
**Status:** `replicated`, `contested`, `capacity-constrained`, `live-tradeable`

---

## Operations

### Ingest a New Source

When the user drops a source in `raw/` (or `../references/`) and says "ingest":

1. **Read** the source in full.
2. **Discuss** key takeaways with the user — what's novel, what contradicts existing wiki pages, what's practically relevant to genesis-alpha.
3. **Write** a paper page in `wiki/papers/` (see format below).
4. **Update** relevant concept pages — add new findings, update empirical magnitudes, flag contradictions with `> ⚠️ Contradicts [[OtherPage]]`.
5. **Update or create** strategy pages if the paper bears on a known strategy approach.
6. **Update** `index.md` — add the new paper page and any new concept/strategy pages.
7. **Append** an ingest entry to `log.md`.

A single paper typically touches 5–15 pages. Prefer updating existing pages to creating new ones unless a concept genuinely lacks its own page.

### Query

When the user asks a question:

1. Read `index.md` to identify relevant pages.
2. Read those pages; synthesize an answer with citations to specific wiki pages.
3. If the answer is substantive (a comparison, multi-paper synthesis, or novel connection), offer to file it as a new page — a good synthesis shouldn't disappear into chat history.

Answers can take forms beyond prose: comparison tables, matplotlib charts (describe the code), Marp slide outlines. Offer these when they'd be more useful.

### Lint

When the user asks for a wiki health check:

1. Read `index.md` and scan all pages.
2. Report: contradictions between pages, stale claims superseded by newer sources, orphan pages (no inbound links), important concepts mentioned but lacking their own page, missing cross-references, data gaps that a web search or new source could fill.
3. Suggest new questions to investigate and new sources to look for.
4. Ask which issues the user wants fixed.
5. Append a lint entry to `log.md`.

---

## Paper Page Format

**Filename:** `wiki/papers/AuthorYYYY-SlugTitle.md` (e.g., `Jegadeesh1993-Momentum.md`)

```markdown
---
type: paper
tags: [signal-type, asset-class]
sources: [AuthorYYYY-SlugTitle]
updated: YYYY-MM-DD
---

# Full Paper Title

**Authors:** Last, F., Last, F.  **Venue / Year:** Journal Name, YYYY

## Contribution
One paragraph — what this paper adds that wasn't established before.

## Methodology
- **Universe:** what instruments, what sample period
- **Signal construction:** exact lookback, formation period, holding period
- **Risk controls:** how they handle beta, sector, size, etc.
- **Benchmark / comparison:** what they compare against

## Key Results
- Quantitative, precise: Sharpe ratios, t-stats, annualized returns, IC
- Note whether figures are pre- or post-cost
- Note any capacity estimates mentioned

## Limitations / Caveats
- Data mining / multiple testing concerns
- Sample period issues (crisis, regime change)
- Capacity, market impact at scale
- Known replication failures or debates

## Connection to genesis-alpha
How does this bear on our alpha/signals/, backtest/costs/, execution/simulator/, or risk/ layer? What would we do differently or test next based on this?

## Related Pages
[[concepts/ConceptName]], [[strategies/StrategyName]], [[papers/RelatedPaper]]
```

---

## Concept Page Format

**Filename:** `wiki/concepts/ConceptName.md` (e.g., `MarketImpact.md`, `ShortTermReversal.md`)

```markdown
---
type: concept
tags: [...]
sources: [slug1, slug2]
updated: YYYY-MM-DD
---

# Concept Name

## Definition
Precise definition. Include units, typical magnitudes, and any important variants.

## Why It Matters for genesis-alpha
Practical relevance to our signals, backtest, execution, or risk layer.

## Empirical Findings
What the literature says, with magnitudes. Call out contradictions explicitly:
> ⚠️ [[papers/Paper1]] finds X; [[papers/Paper2]] finds Y — likely explained by different universes.

## Key Papers
- [[papers/AuthorYYYY]] — one-line summary of that paper's contribution to this concept

## Open Questions
What remains uncertain, contested, or not yet covered in our wiki.

## Related Pages
```

---

## Strategy Page Format

**Filename:** `wiki/strategies/StrategyName.md`

```markdown
---
type: strategy
tags: [signal-type, asset-class]
sources: [slug1, slug2]
updated: YYYY-MM-DD
---

# Strategy Name

## Overview
What this strategy does in plain terms.

## Signal Construction
- Lookback, holding period, weighting scheme
- Long/short or long-only; gross vs net exposure

## Evidence Summary
Consolidated findings across papers: typical Sharpe (pre/post cost), persistence, decay.

## Costs and Capacity
What we know about turnover, market impact, and realistic AUM limits.

## Variants and Extensions
Related signals, enhancements, combinations known to add value.

## genesis-alpha Implementation Notes
Where this would live in the codebase (`alpha/signals/`, `alpha/factors/`), key parameters to tune.

## Related Pages
```

---

## Quant Domain Conventions

- **Sharpe ratios:** always note pre- or post-cost, and whether annualized. Default assumption: annualized pre-cost unless stated.
- **Signal construction:** always note lookback period, holding period, universe, long/short vs long-only.
- **Asset classes:** tag explicitly — the same signal can behave differently across equities vs futures.
- **Replication status:** note if independently replicated, contested, or if it's a single-paper claim.
- **Capacity:** flag `capacity-constrained` tag if paper discusses limits or if signal seems small-cap dependent.
- **Contradictions:** when a new source conflicts with existing pages, don't silently overwrite — add the new finding with a ⚠️ block and note what might explain the divergence (universe, period, methodology).

---

## index.md Format

Organized by page type. Each entry:

```
- [[wiki/papers/AuthorYYYY-Slug|Title]] — one-line summary (N sources)
```

Update on every ingest. The LLM reads this first on every query to find relevant pages.

---

## log.md Format

Each entry header:

```
## [YYYY-MM-DD] operation | description
```

Operations: `ingest`, `query`, `lint`, `update`, `init`

Quick tail: `grep "^## \[" log.md | tail -10`

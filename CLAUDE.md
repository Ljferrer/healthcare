# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Healthcare Marketplace for Claude Code — a collection of 3 skills and 4 MCP server plugins for healthcare workflows. This is not a compiled application; it's a skill/plugin marketplace where the "code" is structured Markdown prompts, JSON configuration, and one Python utility script.

## Repository Structure

**Skills** (each has `SKILL.md` entry point + `references/` subskills + `assets/`):
- `clinical-trial-protocol-skill/` — Multi-step FDA/NIH protocol generator using waypoint architecture
- `prior-auth-review-skill/` — Prior authorization intake + decision workflow (demo)
- `fhir-developer-skill/` — FHIR R4 reference documentation skill (no executable workflow)

**MCP Server Plugins** (each has `.claude-plugin/plugin.json` pointing to remote HTTP endpoint):
- `cms-coverage/` → `mcp.deepsense.ai/cms_coverage/mcp`
- `npi-registry/` → `mcp.deepsense.ai/npi_registry/mcp`
- `icd10-codes/` → `mcp.deepsense.ai/icd10_codes/mcp`
- `pubmed/` → `pubmed.mcp.claude.com/mcp`

**Central config**: `.claude-plugin/marketplace.json` — single source of truth for all plugin registration.

## Key Architecture Patterns

### Waypoint-Based State Persistence
Skills write intermediate results to `waypoints/` as JSON/Markdown files. This allows resuming interrupted workflows and provides an audit trail. Each step reads prior waypoints and writes its own.

### Load-On-Demand Subskills
Skill execution loads subskill files from `references/` one at a time as each step runs. **Do not pre-read all subskill files at startup** — this wastes token budget.

### Human-in-the-Loop
All high-stakes decisions require human confirmation. The prior-auth skill defaults to APPROVE or PEND only (never auto-DENY in lenient mode). Clinical trial protocols are always marked as preliminary drafts requiring professional review.

### Decision Rubric (Prior Auth)
`prior-auth-review-skill/references/rubric.md` controls approval logic. Default is lenient mode (PEND when uncertain). Editable to enable strict mode.

## Commands

**Python dependency** (only needed for clinical-trial-protocol-skill):
```bash
pip install scipy numpy
```

**Run sample size calculator**:
```bash
python clinical-trial-protocol-skill/scripts/sample_size_calculator.py
```

**Install marketplace locally for testing**:
```bash
/plugin marketplace add anthropics/healthcare
```

## CI/CD

GitHub Actions workflows in `.github/workflows/`:
- `claude-skill-review.yml` — On PR, detects changed skills (directories containing `SKILL.md`) and runs AI review
- `claude-code-review.yml` / `claude.yml` — AI review triggered by `@claude` in issue comments

## Adding a New Skill

1. Create `<skill-name>/SKILL.md` with YAML frontmatter (`name`, `description`) and execution instructions
2. Add subskill steps in `<skill-name>/references/`
3. Register in `.claude-plugin/marketplace.json` under `plugins` array with a `skills` path
4. Add sample data in `<skill-name>/assets/` if applicable

## Adding a New MCP Server Plugin

1. Create `<plugin-name>/.claude-plugin/plugin.json` with `mcpServers` defining the HTTP endpoint URL
2. Register in `.claude-plugin/marketplace.json` under `plugins` array (no `skills` field)
3. Document in root `README.md` MCP servers table

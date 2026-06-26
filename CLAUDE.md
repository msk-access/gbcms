# gbcms — Claude Code entry

> Canonical standing context is `AGENTS.md` (cross-tool). This file imports it so
> Claude Code loads the same rules every other tool reads. Keep Claude-only notes
> below the import; put shared rules in `AGENTS.md`.

@AGENTS.md

## Claude-only notes
- Skills load from `.claude/skills/` (a symlink to `.agents/skills/`).
- Hooks: `.claude/settings.json` (time injection + `guard-destructive.sh`).
- Memory auto-recall reads `.agents/memory/` (symlinked from the Claude project memory dir).

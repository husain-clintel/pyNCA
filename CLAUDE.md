# Claude Code Configuration

## Allowed Tools (Autonomous Mode)

Claude may run the following without asking for permission:

### Bash Commands
- `python`, `python3` - Run Python scripts
- `pip`, `pip3` - Install/manage Python packages
- `pytest` - Run tests
- `git` - Version control operations (except force push)
- `ls`, `pwd`, `cd` - Directory navigation
- `mkdir`, `rm`, `cp`, `mv` - File operations
- `cat`, `head`, `tail`, `less` - View files
- `grep`, `find`, `which` - Search utilities
- `npm`, `node` - Node.js operations if needed
- `make` - Build automation

### File Operations
- Read any file in this project
- Edit any file in this project
- Create new files as needed
- Delete files when appropriate

### Restrictions
- Do not push to remote repositories without explicit permission
- Do not modify system files outside this project
- Do not install global packages without permission
- Do not run destructive git commands (force push, hard reset) without permission

## Project Context

This is a Python NCA (Non-Compartmental Analysis) project.

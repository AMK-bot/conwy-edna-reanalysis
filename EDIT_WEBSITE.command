#!/bin/zsh
cd "$(dirname "$0")"
exec /usr/bin/env python3 scripts/edit_website.py

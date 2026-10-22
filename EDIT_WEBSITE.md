# Personalize the website text

This editor changes selected narrative text on the project website without requiring
HTML knowledge. It runs only on your computer and stores changes separately from the
canonical page.

## Open the editor

From the project folder, double-click `EDIT_WEBSITE.command`. If your system does not open
the launcher, run:

```bash
python3 scripts/edit_website.py
```

A browser tab opens automatically. Keep the launcher window open while editing.

## Edit and save

1. Click an outlined heading or paragraph in the live preview.
2. Type or paste plain text.
3. Choose **Save website text**.
4. Use **Open full website** to review the result outside the editor.

**Discard draft** returns to the last saved wording. **Undo last save** restores the
previous saved version. **Restore original** clears all personalizations. **Export copy**
and **Import copy** move the editable wording as a small JSON file.

The editor intentionally locks evidence numbers, accepted/rejected/non-conclusive
decisions, links, figures, and model results. It changes the website only; text embedded
in diagrams and the four-page PDF is not rewritten. Press `Ctrl+C` in the
launcher window to stop the local editor.

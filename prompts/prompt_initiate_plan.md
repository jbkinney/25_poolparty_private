Does this make sense? Any questions, suggestions, concerns, or comments?

Please acknowledge these instructions before proceeding. End by making sure that all tests pass.

Note: I'm using uv to manage this python package, so if you need to run code at the command line do `uv run ...`.

Make the comments and docstrings concise, with listing all indivudal parameters limited to essential user-facing functions.

Finally, write a jupyter notebook that concisely demonstrates this functionality.
- [IGNORE THIS] The notebook should be in notebooks/ and have a name of the form YY.MM.DD_description.ipynb
- Keep the notebook concise
- No markdown cells
- DO NOT ADD COMMENTS; the code and perhaps print statements should speak for themselves.
- If you want the cell to have multiple sections, start each section with a one-line print() statement explaining the section.

When editing Jupyter notebooks (.ipynb files), use the `edit_notebook` tool - do NOT write raw JSON. Key parameters:
- `target_notebook`: path to the notebook
- `cell_idx`: 0-based cell index
- `is_new_cell`: true to create new cell at index, false to edit existing
- `cell_language`: 'python', 'markdown', 'raw', etc.
- `old_string`: text to replace (empty for new cells)
- `new_string`: replacement text or new cell content

To delete cell content, pass empty string as `new_string`. To run notebooks, use `uv run jupyter nbconvert --to notebook --execute <notebook.ipynb> --inplace`.

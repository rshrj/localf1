# Local $\mathbb{F}_1$ Notes
Author: Rishi Raj  
Institution: LPTHE, Sorbonne Université  
Date: October 2025  

## Structure
- `main.tex` — main LaTeX source
- `references/` — BibLaTeX files
- `figures/` — images
- `notebooks/` — Mathematica explorations
- `build/` — compiled artifacts (ignored in Git)

## Build
Run `make build` or

```bash
latexmk -pdf -output-directory=build main.tex
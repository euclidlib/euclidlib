# -- Project information -----------------------------------------------------

project = "euclidlib"
author = "the euclidlib team"

# -- General configuration ---------------------------------------------------

extensions = [
    "myst_nb",               # MyST Markdown + notebooks
    "sphinx.ext.mathjax",    # Math support
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx_design",         # For layout utilities (optional but useful)
]

# MyST NB & MyST settings
myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "tasklist",
    "substitution",
    "linkify",
]

# Notebook execution — no rule specified in myst.yml
# so we keep the default (execute only if needed)
nb_execution_mode = "auto"

# Bibliography (from myst.yml → project.bibliography)
bibtex_bibfiles = ["references.bib"]

# -- HTML theme --------------------------------------------------------------

html_theme = "sphinx_book_theme"

html_logo = "euclidlib-patch.png"     # site.options.logo

html_theme_options = {
    "logo": {
        "image_light": "euclidlib-patch.png",
        "image_dark": "euclidlib-patch.png",
    },
    # "folders": true in myst.yml corresponds to showing folder sidebar structure
    "show_nav_level": 2,
}

# -- PDF Export (Sphinx only partially supports this) ------------------------
# Your myst.yml specified:
#   exports:
#     - format: pdf
#       template: plain_latex_book
#       output: exports/book.pdf
#
# To enable LaTeX builds on ReadTheDocs, you would add:
# latex_engine = "pdflatex"
# latex_elements = {}
#
# But PDF export must be configured in RTD, not here.

###############################################################################
# Cleaned and fixed conf.py for ReadTheDocs + MyST + JupyterBook ecosystem
###############################################################################

author = "the euclidlib team"
copyright = "2024"

extensions = [
    #"myst_parser",
    #"myst_nb",
    "sphinx_togglebutton",
    "sphinx_copybutton",
    "sphinx_thebe",
    "sphinx_comments",
    "sphinx_external_toc",
    "sphinx.ext.intersphinx",
    "sphinx_design",
    "sphinx_book_theme",
    "sphinxcontrib.bibtex",
    "sphinx_jupyterbook_latex",
    "sphinx_multitoc_numbering",
    "sphinx.ext.mathjax",
]

# --------------------------------------------------------------------------
# FIX: Make Sphinx accept Markdown and Notebooks
# --------------------------------------------------------------------------
source_suffix = {
    ".md": "markdown",
    ".rst": "restructuredtext",
    ".ipynb": "myst-nb",
}

exclude_patterns = [
    "**.ipynb_checkpoints",
    ".DS_Store",
    "Thumbs.db",
    "_build",
]

# --------------------------------------------------------------------------
# Jupyter-Book settings
# --------------------------------------------------------------------------
external_toc_path = "_toc.yml"
external_toc_exclude_missing = False

myst_enable_extensions = [
    "colon_fence",
    "dollarmath",
    "linkify",
    "substitution",
    "tasklist",
    "include",
]

myst_url_schemes = ["mailto", "http", "https"]

# --------------------------------------------------------------------------
# Notebook execution – RTD-safe
# (disable or make it cached; RTD often fails with "kernel not found")
# --------------------------------------------------------------------------
nb_execution_mode = "off"
nb_execution_timeout = 30

# --------------------------------------------------------------------------
# HTML theme
# --------------------------------------------------------------------------
html_theme = "sphinx_book_theme"
html_logo = "euclidlib-patch.png"
html_title = "euclidlib"

html_theme_options = {
    "search_bar_text": "Search...",
    "repository_url": "https://github.com/euclidlib/euclidlib",
    "repository_branch": "8-add-readthedocs",
    "path_to_docs": "docs",
    "use_repository_button": True,
    "home_page_in_toc": True,
}

# --------------------------------------------------------------------------
# Bibliography
# --------------------------------------------------------------------------
bibtex_bibfiles = ["references.bib"]

# --------------------------------------------------------------------------
# Latex
# --------------------------------------------------------------------------
latex_engine = "pdflatex"

# --------------------------------------------------------------------------
# Misc
# --------------------------------------------------------------------------
pygments_style = "sphinx"
suppress_warnings = ["myst.domains"]
numfig = True

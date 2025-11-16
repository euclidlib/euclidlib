###############################################################################
# Clean conf.py for ReadTheDocs + JupyterBook 2 + MyST-NB 1.1
###############################################################################

author = "the euclidlib team"
copyright = "2025"

# --------------------------------------------------------------------------
# Extensions — minimal and safe for JB2
# --------------------------------------------------------------------------
extensions = [
    "myst_nb",
    "sphinx_external_toc",
    "sphinx_book_theme",
    "sphinx_design",
    "sphinx_togglebutton",
    "sphinx_copybutton",
    "sphinxcontrib.bibtex",
    "sphinx_thebe",
]

# --------------------------------------------------------------------------
# Never override source_suffix when using MyST-NB (JB2 handles Markdown)
# --------------------------------------------------------------------------
# DO NOT define source_suffix — myst_nb registers markdown automatically.

exclude_patterns = [
    "**.ipynb_checkpoints",
    ".DS_Store",
    "Thumbs.db",
    "_build",
]

# --------------------------------------------------------------------------
# Table of contents from external YAML
# --------------------------------------------------------------------------
external_toc_path = "_toc.yml"
external_toc_exclude_missing = False

# --------------------------------------------------------------------------
# IMPORTANT: Do NOT define myst_enable_extensions or myst_url_schemes.
# Jupyter Book 2 + MyST-Parser 3 configure these automatically.
# --------------------------------------------------------------------------

# --------------------------------------------------------------------------
# Notebook execution (recommended for RTD)
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

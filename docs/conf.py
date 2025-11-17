# Configuration file for the Sphinx documentation builder.

from datetime import date
from sphinx.application import Sphinx

# -- Project information -----------------------------------------------------

project = "euclidlib"
author = "the euclidlib team"
copyright = f"{date.today().year}, {author}"

# Master document
master_doc = "index"
language = "en"

# -- General configuration ---------------------------------------------------

extensions = [
    "myst_nb",  # MyST Markdown + notebooks
    "sphinx.ext.mathjax",  # Math support
    "sphinx.ext.autodoc",  # API documentation
    "sphinx.ext.napoleon",  # NumPy/Google docstrings
    "sphinx_design",  # Layout utilities
]

templates_path = ["_templates"]

exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
]

# -- MyST & Notebook settings ------------------------------------------------

myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "tasklist",
    "substitution",
    "linkify",
]

nb_execution_mode = "auto"  # Execute notebooks only if needed

# -- Bibliography ------------------------------------------------------------

bibtex_bibfiles = ["references.bib"]

# -- HTML output -------------------------------------------------------------

html_theme = "sphinx_book_theme"

html_logo = "euclidlib-patch.png"

html_theme_options = {
    "logo": {
        "image_light": "euclidlib-patch.png",
        "image_dark": "euclidlib-patch.png",
    },
    "show_nav_level": 2,  # Equivalent to "folders: true" in myst.yml
}

html_static_path = ["_static"]
html_css_files = ["local.css"]

# -- PDF / LaTeX output (optional) ------------------------------------------
# If you want ReadTheDocs to build PDF:
# latex_engine = "pdflatex"
# latex_elements = {}

# -- Optional: local Sphinx extensions setup --------------------------------


def setup(app: Sphinx):
    """Custom Sphinx setup for euclidlib."""
    # Add directives, transforms, etc, if needed later.
    pass

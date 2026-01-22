# Configuration file for the Sphinx documentation builder.

from datetime import date
from sphinx.application import Sphinx
import os

# -- Project information -----------------------------------------------------

project = "euclidlib"
author = "the euclidlib team"
copyright = f"{date.today().year}, {author}"

# Master document
language = "en"

# -- General configuration ---------------------------------------------------

extensions = [
    "myst_parser",
    "myst_nb",
    "sphinx.ext.mathjax",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx_design",
]

templates_path = ["_templates"]

exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
]

# Support both RST and Markdown
source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

# -- MyST & Notebook settings ------------------------------------------------

myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "tasklist",
    "substitution",
    "linkify",
]

# Execute notebooks automatically
jupyter_execute_notebooks = "auto"  # RTD executes notebooks if needed

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
    "show_nav_level": 2,  # Depth of sidebar
}

# Only use _static if it exists
if os.path.exists("_static"):
    html_static_path = ["_static"]
    html_css_files = ["local.css"]
else:
    html_static_path = []

# -- PDF / LaTeX output (optional) ------------------------------------------
# latex_engine = "pdflatex"
# latex_elements = {}


# -- Optional: local Sphinx extensions setup --------------------------------
def setup(app: Sphinx):
    """Custom Sphinx setup for euclidlib."""
    pass

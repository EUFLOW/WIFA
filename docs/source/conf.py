# Configuration file for the Sphinx documentation builder.
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

sys.path.insert(0, os.path.abspath("../.."))

# -- Project information -----------------------------------------------------

project = "WIFA"
copyright = "2024, EU-FLOW"
author = "EU-FLOW"

# The short X.Y version
version = "1.0"
# The full version, including alpha/beta/rc tags
release = "1.1"

# GitHub repository info for PyData theme
github_user = "EUFLOW"
github_repo = "WIFA"
github_version = "main"


# -- General configuration ---------------------------------------------------

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.mathjax",
    "sphinx.ext.githubpages",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinxcontrib.bibtex",
    "sphinx_design",
    "sphinx_copybutton",
]

# Napoleon settings for Google/NumPy style docstrings
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = True
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_type_aliases = None

# Autodoc settings
autodoc_default_options = {
    "members": True,
    "member-order": "bysource",
    "special-members": "__init__",
    "undoc-members": True,
    "exclude-members": "__weakref__",
}

# Intersphinx configuration for linking to external documentation
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "xarray": ("https://docs.xarray.dev/en/stable/", None),
    "pandas": ("https://pandas.pydata.org/docs/", None),
}

# Copy button configuration
copybutton_prompt_text = r">>> |\.\.\. |\$ |In \[\d*\]: | {2,5}\.\.\.: | {5,8}: "
copybutton_prompt_is_regexp = True

templates_path = ["_templates"]
bibtex_bibfiles = ["references.bib"]
source_suffix = ".rst"
master_doc = "index"
language = "en"
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

html_theme = "pydata_sphinx_theme"
html_theme_options = {
    "logo": {
        "image_light": "../img/Logo_FLOW.png",
        "image_dark": "../img/Logo_FLOW.png",
        "text": "WIFA",
    },
    "icon_links": [
        {
            "name": "GitHub",
            "url": f"https://github.com/{github_user}/{github_repo}",
            "icon": "fa-brands fa-github",
        },
    ],
    "use_edit_page_button": True,
    "show_toc_level": 2,
    "navigation_with_keys": True,
    "header_links_before_dropdown": 10,
    "show_nav_level": 2,
    "navbar_align": "left",
    "navbar_start": ["navbar-logo"],
    "navbar_center": ["navbar-nav"],
    "navbar_end": ["theme-switcher", "navbar-icon-links"],
    "footer_start": ["copyright"],
    "footer_end": ["last-updated"],
    "pygments_light_style": "default",
    "pygments_dark_style": "monokai",
}

html_context = {
    "github_user": github_user,
    "github_repo": github_repo,
    "github_version": github_version,
    "doc_path": "docs/source",
}

html_static_path = ["_static"]
html_css_files = ["css/custom.css"]
html_favicon = "../img/Logo_FLOW.png"


# -- Options for HTMLHelp output ---------------------------------------------

htmlhelp_basename = "WIFAdoc"

# -- Options for LaTeX output ------------------------------------------------

latex_elements = {}
latex_documents = [
    (master_doc, "WIFA.tex", "WIFA Documentation", "EU-FLOW", "manual"),
]

# -- Options for manual page output ------------------------------------------

man_pages = [(master_doc, "wifa", "WIFA Documentation", [author], 1)]

# -- Options for Texinfo output ----------------------------------------------

texinfo_documents = [
    (
        master_doc,
        "WIFA",
        "WIFA Documentation",
        author,
        "WIFA",
        "Multi-fidelity wind farm simulation framework",
        "Miscellaneous",
    ),
]

# -- Options for Epub output -------------------------------------------------

epub_title = project
epub_exclude_files = ["search.html"]

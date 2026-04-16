import datetime
import inspect
import os
import sys
from contextlib import suppress
from pathlib import Path

# Disable import of swiftest.core so that we don't have to build the Fortran code when building the docs
import swiftest

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
project = "Swiftest"
copyright = f"{datetime.datetime.now().year}, David A. Minton"
author = "David A. Minton"
version = swiftest.__version__
release = version

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.doctest",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "sphinx.ext.extlinks",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "IPython.sphinxext.ipython_directive",
    "IPython.sphinxext.ipython_console_highlighting",
    "nbsphinx",
    "sphinxext.opengraph",
    "sphinx_copybutton",
    "sphinx_design",
    "sphinx_inline_tabs",
]

extlinks = {
    "issue": ("https://github.com/MintonGroup/swiftest/issues/%s", "GH%s"),
    "pull": ("https://github.com/MintonGroup/swiftest/pull/%s", "PR%s"),
    "discussion": ("https://github.com/MintonGroup/swiftest/discussions/%s", "D%s"),
    "release": ("https://github.com/MintonGroup/swiftest/releases/tag/%s", "%s"),
}

# Sometimes the savefig directory doesn't exist and needs to be created
# https://github.com/ipython/ipython/issues/8733
# becomes obsolete when we can pin ipython>=5.2; see ci/requirements/doc.yml
ipython_savefig_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "_build", "html", "_static")
if not os.path.exists(ipython_savefig_dir):
    os.makedirs(ipython_savefig_dir)


templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "**.ipynb_checkpoints"]


# The name of the Pygments (syntax highlighting) style to use.
pygments_style = "sphinx"


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output


# Use Autodoc and Napolean for extracting docstrings
autosummary_generate = True
autodoc_typehints = "none"
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_use_param = False
napoleon_use_rtype = False

autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "show-inheritance": True,
    "member-order": "bysource",
    "no-index-entry": True,
}
autoapi_generate_api_docs = False


intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "xarray": ("https://docs.xarray.dev/en/stable/", None),
}

templates_path = ["_templates"]

html_theme = "sphinx_book_theme"
html_title = ""

html_context = {
    "github_user": "profminton",
    "github_repo": "swiftest",
    "github_version": "main",
    "doc_path": "docs",
}

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
html_theme_options = {
    # analytics_id=''  this is configured in rtfd.io
    # canonical_url="",
    "repository_url": "https://github.com/MintonGroup/swiftest",
    "repository_branch": "master",
    "navigation_with_keys": False,  # pydata/pydata-sphinx-theme#1492
    "path_to_docs": "docs",
    "use_edit_page_button": True,
    "use_repository_button": True,
    "use_issues_button": True,
    "home_page_in_toc": False,
    "extra_footer": """Theme by the <a href="https://ebp.jupyterbook.org">Executable Book Project</a></p>""",
    "icon_links": [],  # workaround for pydata/pydata-sphinx-theme#1220
    "announcement": "🍾 <a href='https://github.com/MintonGroup/swiftest/discussions/1'>Swiftest is currently under development</a> 🎉",
    "show_toc_level": 4,
    "collapse_navigation": True,
    "navigation_depth": 2,
}


# The name of an image file (relative to this directory) to place at the top
# of the sidebar.
html_logo = "_static/logos/swiftest_social_preview.svg"

# The name of an image file (within the static path) to use as favicon of the
# docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
html_favicon = "_static/logos/swiftest_logo.svg"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]
html_css_files = ["style.css"]

# configuration for sphinxext.opengraph
ogp_site_url = "https://swiftest.readthedocs.io/en/latest/"
ogp_image = "https://swiftest.readthedocs.io/en/stable/_static/logos/swiftest_social_preview.png"
ogp_custom_meta_tags = [
    '<meta name="image" property="og:image" content="https://swiftest.readthedocs.io/en/stable/_static/logos/swiftest_social_preview.png" />',
]


def html_page_context(app, pagename, templatename, context, doctree):
    # Disable edit button for docstring generated pages
    if "generated" in pagename:
        context["theme_use_edit_page_button"] = False

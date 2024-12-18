# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'gebt'
copyright = '2024, Wenbin Yu'
author = 'Wenbin Yu'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.mathjax',
    'sphinx.ext.napoleon',
    "sphinx.ext.githubpages",
    'sphinx-prompt',
    'sphinx_copybutton',
]

templates_path = ['_templates']
exclude_patterns = []

html_theme_options = {
    'logo': {
        'text': 'GEBT',
    },
    'show_nav_level': 2,
    # "path_to_docs": "doc/source",
    # 'use_edit_page_button': True,
    # "use_repository_button": True,
    # "use_issues_button": True,
    # 'collapse_navigation': True,
    'navigation_depth': 4,
    "announcement": "Documentation is under construction.",
}



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'pydata_sphinx_theme'
html_static_path = ['_static']

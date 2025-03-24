# -*- coding: utf-8 -*-

import sys

extensions = [
    # 'sphinx.ext.autodoc',
    'sphinx.ext.mathjax',
    'sphinx.ext.intersphinx',
    'matplotlib.sphinxext.plot_directive',
    # 'sphinx.ext.doctest',
    'sphinx.ext.todo',
    # 'sphinx.ext.viewcode',
    # 'sphinx.ext.autosummary',
]

source_suffix = '.rst'
todo_include_todos = True

project = u'ChiQ'
copyright = u'2025-, The University of Tokyo'
author = u'ChiQ Developers'

html_theme = 'wild'
import wild_sphinx_theme
html_theme_path = [wild_sphinx_theme.get_theme_dir()]

html_show_sphinx = False
html_context = {'header_title': 'ChiQ',
                'header_subtitle': 'χ(q) solver based on the Bethe-Salpeter equation',
                'header_links': [['Install', 'install'],
                                 ['Documentation', 'documentation'],
                                 #['Presentatation', 'presentation'],
                                 ['Issues', 'issues'],
                                 ['About ChiQ', 'about']]}
#html_static_path = ['@CMAKE_SOURCE_DIR@/doc/_static']
html_static_path = ['_static']
#html_sidebars = {'index': ['sideb.html', 'searchbox.html']}
html_sidebars = {'**': ['globaltoc.html', 'relations.html', 'searchbox.html']}

#html_title = "ChiQ software"
#html_logo = 'logo_chiq.png'

# no 'module' in the header
html_domain_indices = False

# no 'index' in the header
html_use_index = False

htmlhelp_basename = 'ChiQdoc'

rst_epilog = '.. |BRANCH| replace:: {}'.format('@GIT_BRANCH_NAME@')

latex_engine = "lualatex"

# overwrite css of html_theme
def setup(app):
    app.add_css_file('chiq.css')

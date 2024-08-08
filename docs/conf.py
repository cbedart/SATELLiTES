from pathlib import Path
CURR_PATH = Path(__file__).absolute().parent

project = 'SATELLiTES'
copyright = '2024, Corentin BEDART'
author = 'Corentin BEDART'
release = '1.0.9'

extensions = []
templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

html_theme = 'sphinx_rtd_theme'
html_theme_options = {
    'includehidden': False,
    'logo_only': True,
}

html_static_path = ['_static']
html_logo = str(CURR_PATH / 'SATELLiTES_favicon.png')

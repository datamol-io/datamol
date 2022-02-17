# Configuration

$PROJECT = $GITHUB_REPO  = 'datamol'
$GITHUB_ORG = 'datamol-org'
$PUSH_TAG_REMOTE = 'git@github.com:datamol-org/datamol.git'
$GHRELEASE_TARGET = 'main'

# Logic

$AUTHORS_FILENAME = 'AUTHORS.rst'
$AUTHORS_METADATA = '.authors.yml'
$AUTHORS_SORTBY = 'alpha'
$AUTHORS_MAILMAP = '.mailmap'

$CHANGELOG_FILENAME = 'CHANGELOG.rst'
$CHANGELOG_TEMPLATE = 'TEMPLATE.rst'
$CHANGELOG_NEWS = 'news'

$PYPI_BUILD_COMMANDS = ['sdist']
$PYPI_UPLOAD = True

# Current pypi action is broken. See https://github.com/regro/rever/issues/254
# For now you must upload manually: `python setup.py sdist && twine upload dist/*`
# $ACTIVITIES = ['check', 'authors', 'changelog', 'version_bump', 'tag', 'push_tag', 'ghrelease', 'pypi']
$ACTIVITIES = ['check', 'authors', 'changelog', 'version_bump', 'tag', 'push_tag', 'ghrelease']

$VERSION_BUMP_PATTERNS = [('datamol/_version.py', r'__version__\s*=.*', "__version__ = \"$VERSION\""),
                          ('setup.py', r'version\s*=.*,', "version=\"$VERSION\",")
                          ]

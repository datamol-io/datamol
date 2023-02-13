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

$ACTIVITIES = ['check', 'authors', 'changelog', 'version_bump', 'tag', 'push_tag', 'ghrelease']

$VERSION_BUMP_PATTERNS = [('datamol/_version.py', r'__version__\s*=.*', "__version__ = \"$VERSION\""),
                          ('pyproject.toml', r'version\s*=.*  # project', "version = \"$VERSION\"  # project")
                          ]

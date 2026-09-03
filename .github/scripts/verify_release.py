"""Validate manual release inputs before tests or publication."""

import os
import re
from pathlib import Path

VERSION = re.compile(r"(0|[1-9]\d*)\.(0|[1-9]\d*)\.(0|[1-9]\d*)(?:(a|b|rc)(0|[1-9]\d*))?")


def validate_release(version, ref, dry_run, changelog):
    """Return whether a canonical release version is a prerelease."""
    match = VERSION.fullmatch(version)
    if match is None:
        raise ValueError("Use a version such as 1.0.0 or 1.0.0rc1, without a v prefix.")
    if not dry_run:
        if ref != "refs/heads/main":
            raise ValueError("Publication must be launched from main.")
        heading = rf"^## \[?{re.escape(version)}\]?(?:[ \t].*)?$"
        if not re.search(heading, changelog, flags=re.MULTILINE):
            raise ValueError(
                f"Finalize CHANGELOG.md with a '## {version} - YYYY-MM-DD' heading first."
            )
    return match.group(4) is not None


if __name__ == "__main__":
    version = os.environ["RELEASE_VERSION"]
    prerelease = validate_release(
        version,
        os.environ["GITHUB_REF"],
        os.environ["DRY_RUN"] == "true",
        Path("CHANGELOG.md").read_text(),
    )
    with open(os.environ["GITHUB_OUTPUT"], "a") as output:
        output.write(f"version={version}\nprerelease={str(prerelease).lower()}\n")

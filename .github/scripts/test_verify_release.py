import unittest

from verify_release import validate_release


class ReleaseValidationTests(unittest.TestCase):
    def test_stable_and_prerelease_versions(self):
        for version in ("1.0.0", "3.0.0", "1.2.3rc1", "1.2.3a0", "1.2.3b2"):
            with self.subTest(version=version):
                result = validate_release(
                    version, "refs/heads/main", False, f"## {version} - 2026-09-02"
                )
                self.assertEqual(result, any(char in version for char in "abc"))

    def test_rejects_invalid_or_ambiguous_versions(self):
        for version in (
            "",
            "v1.0.0",
            "1.0",
            "01.0.0",
            "1.0.0.dev1",
            "1.0.0+local",
            "1.0.0\n",
            "1.0.0; exit",
        ):
            with self.subTest(version=version), self.assertRaises(ValueError):
                validate_release(version, "refs/heads/main", True, "")

    def test_publication_requires_main(self):
        for ref in ("refs/heads/dev", "refs/tags/1.0.0", "refs/heads/main-other"):
            with self.subTest(ref=ref), self.assertRaisesRegex(ValueError, "main"):
                validate_release("1.0.0", ref, False, "## 1.0.0")

    def test_publication_requires_matching_release_notes(self):
        for notes in ("## Next major release (unreleased)", "## 1.0.01", "## 1.0.0rc1"):
            with self.subTest(notes=notes), self.assertRaisesRegex(ValueError, "CHANGELOG"):
                validate_release("1.0.0", "refs/heads/main", False, notes)

    def test_dry_run_allows_dev_and_unreleased_notes(self):
        self.assertFalse(validate_release("1.0.0", "refs/heads/dev", True, "unreleased"))

    def test_linked_release_heading(self):
        self.assertFalse(
            validate_release("1.0.0", "refs/heads/main", False, "## [1.0.0] - 2026-09-02")
        )


if __name__ == "__main__":
    unittest.main()

"""Keep topic-organized documentation reachable and included in NSE exports."""
from pathlib import Path
import re
import sys
import unittest

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from yaml_support import load_yaml


class DocumentLayoutTests(unittest.TestCase):
    def test_document_roots_only_contain_index(self):
        for library in ("NSE", "ReactingFlow"):
            directory = ROOT / "SolverLibrary" / library / "docs"
            self.assertEqual(sorted(p.name for p in directory.glob("*.md")), ["README.md"])

    def test_nse_manifests_include_topic_documents(self):
        directory = ROOT / "SolverLibrary/NSE"
        expected = {p.relative_to(directory).as_posix() for p in (directory / "docs").rglob("*.md")}
        for filename, component in (("solver_manifest.yaml", "runtime_tools"),
                                    ("solver_manifest_multicomponent.yaml", "foundation_docs")):
            manifest = load_yaml(directory / filename)
            actual = set(manifest["components"][component]["files"])
            self.assertTrue(expected <= actual, (filename, sorted(expected - actual)))
            for path in actual:
                if path.startswith("docs/"):
                    self.assertTrue((directory / path).is_file(), (filename, path))

    def test_local_markdown_links(self):
        paths = [ROOT / "README.md", ROOT / "SolverLibrary/README.md"]
        for directory in (ROOT / "docs", ROOT / "SolverLibrary/NSE/docs",
                          ROOT / "SolverLibrary/ReactingFlow/docs"):
            paths.extend(directory.rglob("*.md"))
        for path in paths:
            for target in re.findall(r"\]\(([^\s)]+)\)", path.read_text(encoding="utf-8")):
                if re.match(r"^(?:[a-zA-Z]+:|/|#)", target):
                    continue
                destination = target.split("#", 1)[0]
                self.assertTrue((path.parent / destination).exists(), (str(path), target))


if __name__ == "__main__":
    unittest.main()

"""Create a zip archive of the repository excluding gitignored files."""

import argparse
import subprocess
import zipfile
from pathlib import Path


def list_repo_files() -> list[Path]:
    """Return tracked and untracked (non-ignored) files using git."""
    result = subprocess.run(
        ["git", "ls-files", "-z", "--cached", "--others", "--exclude-standard"],
        check=True,
        capture_output=True,
        text=True,
    )
    files = [Path(p) for p in result.stdout.split("\0") if p]
    gitignore = Path(".gitignore")
    if gitignore.exists() and gitignore not in files:
        files.append(gitignore)
    return files


def latest_commit_date() -> str:
    """Latest commit date in YYYY-MM-DD, or 'no-commit' if unavailable."""
    try:
        result = subprocess.run(
            ["git", "log", "-1", "--format=%cs"],
            check=True,
            capture_output=True,
            text=True,
        )
        date = result.stdout.strip()
        return date or "no-commit"
    except subprocess.CalledProcessError:
        return "no-commit"


def make_zip(paths: list[Path], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(output, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for path in paths:
            zf.write(path, arcname=path.as_posix())


def main() -> None:
    commit_date = latest_commit_date()
    repo_name = Path(".").resolve().name

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path(f"tmp/snapshot_{repo_name}_{commit_date}.zip"),
        help="Output zipfile path (default: tmp/snapshot_<repo>_<commit-date>.zip)",
    )
    args = parser.parse_args()

    repo_files = list_repo_files()
    make_zip(repo_files, args.output)
    print(f"Archive written to {args.output}")


if __name__ == "__main__":
    main()

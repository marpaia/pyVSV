# Releasing pyVSV

The PyPI distribution is `pyVSV`; the Python import and command are `vsv`. The version is defined in `pyproject.toml`, and the Python API and CLI read the installed distribution metadata.

## Distribution contents

The wheel includes the package, typing marker, ontology, and reference SBOL file. The curated design library remains a separate directory supplied to `vsv.library(path)`. The source archive additionally includes tests, documentation, reference data, and curated designs so its test suite can run independently of a checkout.

## Validate a release

Version changes update `pyproject.toml` and `uv.lock`, with release notes recorded in `CHANGELOG.md`. Validation runs from the repository root:

```sh
uv lock
make setup
make check
version=$(uv run --locked python -c 'import vsv; print(vsv.__version__)')
make check-dist DIST_DIR="dist/$version"
```

Use an empty output directory for the version being prepared so older artifacts are not included. `uv build` builds the wheel from the source distribution. `make check-dist` also checks metadata and README rendering with strict Twine validation.

Install each artifact in a fresh environment and check it without importing from the checkout:

```sh
validation_root=$(mktemp -d)
for artifact in "dist/$version"/*.whl "dist/$version"/*.tar.gz; do
  environment="$validation_root/$(basename "$artifact")"
  uv venv "$environment" --python 3.12
  uv pip install --python "$environment/bin/python" "$artifact"
  "$environment/bin/python" -I scripts/check_install.py
done
```

The install check verifies version agreement, license metadata, the typing marker, packaged resources, public subpackages, and both CLI entry points. The source archive also includes everything needed to run the existing tests:

```sh
source_dir="$validation_root/source"
mkdir "$source_dir"
tar -xzf "dist/$version"/*.tar.gz -C "$source_dir" --strip-components=1
uv run --directory "$source_dir" --locked pytest
```

CI runs the checks on Python 3.12, 3.13, and 3.14 before building distributions. The build job checks both isolated installs and runs the source archive's tests on Python 3.12. It uploads the validated wheel and source archive as the `distributions` artifact.

## Publish

1. Commit and push the release changes, then confirm CI passes for that exact commit.
2. Create and push an annotated tag matching the package version exactly, such as `v0.3.0`. Verify the tag points to the intended commit.
3. Create a GitHub release at that tag, using the corresponding section of `CHANGELOG.md` as release notes. Draft releases do not upload packages.
4. Publish the GitHub release. `Publish to PyPI` runs the same CI workflow again, rejects a tag that differs from the package version, and publishes the validated `distributions` artifact through the `pypi` environment. Approve the environment deployment if required reviewers are configured.
5. Confirm both files appear on [PyPI](https://pypi.org/project/pyVSV/) and install that exact version in a fresh environment:

   ```sh
   uv venv "$validation_root/pypi" --python 3.12
   uv pip install --python "$validation_root/pypi/bin/python" --index-url https://pypi.org/simple "pyVSV==$version"
   "$validation_root/pypi/bin/python" -I scripts/check_install.py
   ```

Publishing a GitHub release triggers a real PyPI upload, including releases marked as prereleases. Use a PEP 440 prerelease version such as `0.3.0rc1` with a matching tag when intended. PyPI does not allow replacing uploaded files; corrections require a new version.

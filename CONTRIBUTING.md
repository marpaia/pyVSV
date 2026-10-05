# Development

Install [uv](https://docs.astral.sh/uv/getting-started/installation/), then run `make setup` to create the local environment from `uv.lock`. The default interpreter is Python 3.12, pinned in `.python-version`.

## Working on the code

Genomes live in `vsv.genome`. Architectures, edits, evidence, and the design library live in `vsv.design`. The command-line interface lives in `src/vsv/cli/`. A default `uv sync` installs the library and the dev tools. Add a dependency when a feature needs it.

Use type annotations for public interfaces. Keep file access, network access, and command-line argument handling separate from analysis code, so computations can be tested independently. Imports should not initiate downloads, write files, or run analyses.

Put automated tests in `tests/`. Use small synthetic fixtures for sequence-related tests, and document their purpose. Tests must run offline without external services or large datasets. Test observable behavior, including error cases, when introducing functionality.

Use notebooks for exploration and move reusable logic into the package. Clear notebook outputs before committing and keep generated artifacts in `results/`.

## Dependencies

```sh
uv add <package>
uv add --dev <development-tool>
uv lock --upgrade-package <package>
```

Review and include both `pyproject.toml` and `uv.lock` whenever changing dependencies. Normal development and CI use `--locked` so stale dependency declarations fail instead of silently changing the lockfile.

## Checks

```sh
make format       # Apply Ruff formatting
make check        # Lint, formatting check, strict type checking, and tests
make build        # Build a source distribution and wheel in dist/
```

Individual checks are available through `make lint`, `make format-check`, `make typecheck`, and `make test`. These targets run tools through uv, so activating the virtual environment is optional.

CI runs these checks on Python 3.12, 3.13, and 3.14, then installs the built wheel into a separate environment to exercise its entry points. The distribution name is `vsv-bio` and its version is defined in `pyproject.toml`. The Python API and CLI read that installed distribution metadata. The import package and the `vsv` command remain `vsv`.


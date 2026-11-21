# Contributing to arcwind

Thank you for your interest in arcwind! Bug reports, questions, documentation
fixes, new examples and code are all welcome. This file explains how to report
a problem, set up a development copy, and propose a change. Support and
decision-making are described in [SUPPORT.md](SUPPORT.md) and
[GOVERNANCE.md](GOVERNANCE.md). Everyone taking part is expected to follow the
[Code of Conduct](CODE_OF_CONDUCT.md).

## Asking questions and reporting problems

Please use the [GitHub issue tracker](https://github.com/claw-astro/arcwind/issues)
for bugs, feature requests and questions about using the package. Search the
open and closed issues first: your question may already have an answer.

A good bug report includes:

- the arcwind version (`python -c "import arcwind; print(arcwind.__version__)"`),
  your Python version and operating system;
- what you ran (the configuration file and the few lines of Python or the
  `arcwind` command), what you expected, and what happened, with the full error
  message;
- if you can, a reproduction on one of the synthetic cubes
  (`arcwind.mock.MockWind.elias20_like()` or `.hh30_like()`), so the problem can
  be investigated without your data. Please don't attach large data cubes to an
  issue; describe the cube (size, channel width, beam) instead.

## Setting up a development copy

```bash
git clone https://github.com/claw-astro/arcwind.git
cd arcwind
pip install -e ".[all,test,docs]"
pytest -q                      # the test suite should pass before you start
```

Working in a fresh virtual or conda environment is recommended.

## Proposing a change

1. For anything beyond a small fix, open an issue first to discuss the idea. This
   avoids duplicated work and lets us agree on the approach (for example, the
   form of a new velocity law or a new configuration option).
2. Fork the repository and create a branch from `main`
   (`git checkout -b fix-noise-estimate`).
3. Make your change, with tests and documentation (see below). Keep each pull
   request focused on one change.
4. Run the checks locally:
   ```bash
   pytest -q
   ruff check --select F,E9 src tests examples
   sphinx-build -E -a -b html docs docs/_build/html   # if you changed the docs
   ```
5. Add a line describing the change under **Unreleased** at the top of
   `docs/changelog.rst`.
6. Open a pull request against `main` and describe what changed and why, with a
   link to the issue it addresses. Continuous integration runs the tests on
   every pull request. The maintainer will review it, may ask for changes, and
   merges it once it is ready.

## Guidelines

**Tests.** New features and bug fixes come with tests in `tests/`. Tests should
be fast and self-contained: build inputs with `arcwind.mock` rather than reading
real data, and keep random seeds fixed so results are reproducible.

**Code.**
- Follow the existing style (PEP 8, line length about 95 characters) and keep
  the package working on Python 3.10, the oldest version supported.
- Put units in names where they matter (`z_au`, `vz_kms`, `y0_arcsec`).
- New options belong in the configuration classes in `config.py`, each with a
  default, a unit and a one-line description; the configuration reference in the
  docs is generated from them.
- New velocity laws are added with `@register_law` (see the "custom law"
  tutorial) rather than by editing the fitter.

**Docstrings and documentation.**
- Use NumPy-style docstrings for every public class and function.
- Write in American spelling.
- In reStructuredText, put expressions with vertical bars in double backquotes
  (``` ``|j|`` ```); a bare `|j|` is read as a substitution.
- Figures used in the docs live in `docs/_static/` and are made by the scripts in
  `examples/`; regenerate them with the script rather than editing the images.

**Scientific changes.** Changes that alter fit results (a new default, a
different noise estimate, a change to the model) should say so in the pull
request and the changelog, with a before/after comparison on a mock cube.

## AI-assisted contributions

Contributions written with the help of AI tools are welcome, on the same terms
as any other contribution. Say in the pull request which parts were AI-assisted.
You are responsible for reviewing, testing and understanding everything you
submit, and you must have the right to contribute it under the MIT license.

## Releases

Releases are made by the maintainer. Versions follow
[semantic versioning](https://semver.org): `MAJOR.MINOR.PATCH`, with breaking
changes only in a new minor version while the package is below 1.0, and
announced in the changelog. Each release is a git tag (`vX.Y.Z`) with a
changelog entry.

## License

By contributing, you agree that your contributions are released under the
project's [MIT License](LICENSE).

# Governance

## Roles

arcwind is currently a single-maintainer project.

- **Maintainer:** Charles J. Law (University of Minnesota). The maintainer reviews
  and merges pull requests, triages issues, decides the scope and direction of
  the package, makes releases, and enforces the [Code of Conduct](CODE_OF_CONDUCT.md).
- **Contributors:** anyone who opens an issue, comments on a discussion,
  improves the documentation or submits a pull request. Contributors are
  credited in the changelog and in the release notes.

## How decisions are made

- Day-to-day decisions (bug fixes, documentation, small features) are made by the
  maintainer through the normal review of issues and pull requests.
- Larger changes, such as a new model component, a change to defaults that
  alters fit results, or a change to the public interface, are discussed first in
  a public GitHub issue. Anyone may comment. The maintainer makes the final
  decision and records the reasoning in the issue, so decisions and their
  rationale stay public.
- Scientific choices (the model, its parameters and the fitting method) are
  documented in the docs and the changelog, so results can be traced to the
  version that produced them.

## Becoming a maintainer

Contributors who have made sustained, high-quality contributions may be invited
to become maintainers. With more than one maintainer, pull requests by a
maintainer are reviewed by another maintainer before merging, and decisions
that cannot be agreed in discussion are made by the original maintainer. This
document will be updated when the project gains maintainers.

## Releases and compatibility

- Releases follow [semantic versioning](https://semver.org). While the package
  is below 1.0, breaking changes are allowed in a new minor version
  (0.x → 0.(x+1)) and are always listed in the changelog.
- Every release is a tagged commit (`vX.Y.Z`) with a changelog entry.
- Removed or renamed features are announced in the changelog.

## Continuity

The code is released under the MIT License, so anyone may continue or fork the
project. If the maintainer can no longer maintain arcwind, the repository will be
transferred to a new maintainer or to a group organization, or clearly marked as
unmaintained, so that users know its status.

## Changes to this document

Changes to governance are proposed in a pull request and announced in the
changelog.

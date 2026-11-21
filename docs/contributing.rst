Contributing, support and governance
====================================

arcwind is developed in the open on
`GitHub <https://github.com/claw-astro/arcwind>`_. Bug reports, questions,
documentation fixes and code are all welcome.

Getting help
------------

Read the documentation first, in particular the :doc:`faq`. For questions,
bugs and feature requests, open an issue on the
`issue tracker <https://github.com/claw-astro/arcwind/issues>`_. Please use the
tracker rather than private email, so that answers are public and help other
users. If you can, reproduce the problem on a synthetic cube
(:class:`arcwind.mock.MockWind`), so it can be investigated without your data.

Support is provided on a best-effort basis by the maintainer, who aims to
respond to new issues within two weeks. Bugs that give wrong results take
priority over new features, and fixes go into the latest release. The full
statement is in
`SUPPORT.md <https://github.com/claw-astro/arcwind/blob/main/SUPPORT.md>`_.

Contributing
------------

.. code-block:: bash

   git clone https://github.com/claw-astro/arcwind.git
   cd arcwind
   pip install -e ".[all,test,docs]"
   pytest -q

Open an issue to discuss anything beyond a small fix, then send a pull request
from a branch, with tests, documentation and a changelog line. Continuous
integration runs the tests on every pull request. The guidelines (code style,
tests, docstrings, scientific changes and AI-assisted contributions) are in
`CONTRIBUTING.md <https://github.com/claw-astro/arcwind/blob/main/CONTRIBUTING.md>`_.

Governance
----------

arcwind is currently maintained by Charles J. Law (University of Minnesota),
who reviews and merges contributions, makes releases and decides the direction
of the package. Larger changes, such as new model components, changes to
defaults that alter fit results or to the public interface, are discussed in a
public issue before they are made, and the decision is recorded there. Releases
follow semantic versioning, each with a git tag and a changelog entry. See
`GOVERNANCE.md <https://github.com/claw-astro/arcwind/blob/main/GOVERNANCE.md>`_
and the
`Code of Conduct <https://github.com/claw-astro/arcwind/blob/main/CODE_OF_CONDUCT.md>`_.

Acknowledgments
---------------

Portions of this package—including initial scaffolding, unit tests, code
refactoring, figure generation, and documentation—were developed with
assistance from Claude (Opus 5.5) by Anthropic. All code has been reviewed and
tested by the maintainers. The human authors reviewed, edited and validated all
AI-assisted outputs, and made the core design decisions: the scientific method,
the wind model and its parameterization, the fitting and uncertainty strategy,
and the structure and interfaces of the package.

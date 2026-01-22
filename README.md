# euclidlib

[![PyPI version](https://img.shields.io/pypi/v/euclidlib.svg)](https://pypi.org/project/euclidlib/)
[![CI](https://github.com/euclidlib/euclidlib/actions/workflows/tests.yml/badge.svg?branch=main)](https://github.com/euclidlib/euclidlib/actions/workflows/tests.yml)
[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white)](https://pre-commit.com/)
[![Tests: pytest](https://img.shields.io/badge/tests-pytest-blue?logo=pytest)](https://docs.pytest.org/)
[![Linting: Ruff](https://img.shields.io/badge/linting-ruff-purple?logo=ruff&logoColor=white)](https://docs.astral.sh/ruff/)
[![Code Style: Prettier](https://img.shields.io/badge/code%20style-prettier-ff69b4.svg?logo=prettier&logoColor=white)](https://prettier.io/)
[![Type Checking: mypy](https://img.shields.io/badge/type%20checking-mypy-8A2BE2?logo=mypy&logoColor=white)](https://mypy.readthedocs.io/)
[![All Contributors](https://img.shields.io/github/all-contributors/euclidlib/euclidlib?color=ee8449&style=flat-square)](#contributors)

## Table of Contents

- [Introduction](#introduction)
- [Installation](#installation)
- [Contributing](#contributing)
- [License](#license)
- [Contributors](#contributors)

## Introduction

`euclidlib` is an unofficial Python package designed to access official Euclid mission products provided by the Science Ground Segment. Its goal is to offer the Euclid community a user-friendly, ready-to-use library that enables immediate work with science-ready Euclid data.

The package is maintained on a best-effort basis by volunteers and contributors within the Euclid community. See the contributor list below.

## Installation

As simple as:

```sh
pip install euclidlib
```

### Prerequisites

- `python>3.7`
- `fitsio`
- `numpy`

## Structure and Format of `euclidlib`

The design of the `euclidlib` package closely follows the organisation of the [Euclid Data Product Description Documentation](http://st-dm.pages.euclid-sgs.uk/data-product-doc/dm10/) and reflects the structure of the Euclid Science Ground Segment.

```mermaid
graph TD
    EUCLIDLIB[euclidlib]

    LE3[le3]
    PK[pk_wl]
    TWOPCF[twopcf_wl]

    PHZ[phz]

    EUCLIDLIB --> LE3
    EUCLIDLIB --> PHZ

    LE3 --> PK
    LE3 --> TWOPCF
```

`euclidlib` provides all data products in a unified, Pythonic format based on dataclasses, ensuring consistent, intuitive, and easy-to-use interfaces across all supported products. Please consult the full documentation for additional details.

## Contributing

If you would like to contribute, follow the following steps:

1. Open an issue to let the `euclidlib` maintainers know about your contribution plans (new Euclid product? New feature? A suggestion?)
2. Create a new branch:
   ```sh
   git checkout -b feature/your-feature-name
   ```
3. Commit your changes:
   ```sh
   git commit -m 'Add some feature'
   ```
4. Push to the branch:
   ```sh
   git push origin feature/your-feature-name
   ```
5. Open a pull request

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contributors

This project follows the [all-contributors](https://github.com/all-contributors/all-contributors) specification. Contributions of any kind are welcome!

To discover the meaning of each icon, hover your mouse over it.

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->
<table>
  <tbody>
    <tr>
      <td align="center" valign="top" width="14.28%"><a href="http://gcanasherrera.com"><img src="https://avatars.githubusercontent.com/u/13239454?v=4?s=100" width="100px;" alt="Guadalupe Cañas-Herrera"/><br /><sub><b>Guadalupe Cañas-Herrera</b></sub></a><br /><a href="#code-gcanasherrera" title="Code">💻</a> <a href="#review-gcanasherrera" title="Reviewed Pull Requests">👀</a> <a href="#ideas-gcanasherrera" title="Ideas, Planning, & Feedback">🤔</a> <a href="#maintenance-gcanasherrera" title="Maintenance">🚧</a> <a href="#test-gcanasherrera" title="Tests">⚠️</a> <a href="#example-gcanasherrera" title="Examples">💡</a></td>
      <td align="center" valign="top" width="14.28%"><a href="http://ntessore.page"><img src="https://avatars.githubusercontent.com/u/3993688?v=4?s=100" width="100px;" alt="Nicolas Tessore"/><br /><sub><b>Nicolas Tessore</b></sub></a><br /><a href="#code-ntessore" title="Code">💻</a> <a href="#review-ntessore" title="Reviewed Pull Requests">👀</a> <a href="#ideas-ntessore" title="Ideas, Planning, & Feedback">🤔</a> <a href="#example-ntessore" title="Examples">💡</a> <a href="#maintenance-ntessore" title="Maintenance">🚧</a> <a href="#test-ntessore" title="Tests">⚠️</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/zahrabaghkhani"><img src="https://avatars.githubusercontent.com/u/47903409?v=4?s=100" width="100px;" alt="Zahra Baghkhani"/><br /><sub><b>Zahra Baghkhani</b></sub></a><br /><a href="#code-zahrabaghkhani" title="Code">💻</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://jaimeruizzapatero.net/"><img src="https://avatars.githubusercontent.com/u/39957598?v=4?s=100" width="100px;" alt="Jaime RZ"/><br /><sub><b>Jaime RZ</b></sub></a><br /><a href="#review-JaimeRZP" title="Reviewed Pull Requests">👀</a> <a href="#ideas-JaimeRZP" title="Ideas, Planning, & Feedback">🤔</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/itutusaus"><img src="https://avatars.githubusercontent.com/u/20775836?v=4?s=100" width="100px;" alt="itutusaus"/><br /><sub><b>itutusaus</b></sub></a><br /><a href="#review-itutusaus" title="Reviewed Pull Requests">👀</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/FelicitasKeil"><img src="https://avatars.githubusercontent.com/u/70713596?v=4?s=100" width="100px;" alt="Felicitas Keil"/><br /><sub><b>Felicitas Keil</b></sub></a><br /><a href="#code-FelicitasKeil" title="Code">💻</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/WillHartley"><img src="https://avatars.githubusercontent.com/u/6814229?v=4?s=100" width="100px;" alt="WillHartley"/><br /><sub><b>WillHartley</b></sub></a><br /><a href="#ideas-WillHartley" title="Ideas, Planning, & Feedback">🤔</a> <a href="#data-WillHartley" title="Data">🔣</a></td>
    </tr>
    <tr>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/FlorianDubath"><img src="https://avatars.githubusercontent.com/u/9742907?v=4?s=100" width="100px;" alt="FlorianDubath"/><br /><sub><b>FlorianDubath</b></sub></a><br /><a href="#ideas-FlorianDubath" title="Ideas, Planning, & Feedback">🤔</a> <a href="#data-FlorianDubath" title="Data">🔣</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/jacopo-salvalaggio"><img src="https://avatars.githubusercontent.com/u/99494103?v=4?s=100" width="100px;" alt="Jacopo Salvalaggio"/><br /><sub><b>Jacopo Salvalaggio</b></sub></a><br /><a href="#code-jacopo-salvalaggio" title="Code">💻</a> <a href="#ideas-jacopo-salvalaggio" title="Ideas, Planning, & Feedback">🤔</a> <a href="#data-jacopo-salvalaggio" title="Data">🔣</a></td>
    </tr>
  </tbody>
</table>

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->

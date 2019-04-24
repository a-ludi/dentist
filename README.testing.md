Testing DENTIST
===============

[![standard-readme compliant](https://img.shields.io/badge/readme%20style-standard-brightgreen.svg?style=flat-square)](https://github.com/RichardLitt/standard-readme)

> Create ground-truth scenarios and check results of gap closing software.

In order to assess the performance of gap closing software, a ground-truth scenario is generated where the true sequence in gaps is know. The results can
be evaluated automatically.


Table of Contents
-----------------

- [Install](#install)
- [Usage](#usage)
- [Maintainer](#maintainer)
- [Contributing](#contributing)
- [License](#license)


Install
--------

Install DENTIST as usual but pass the additional flag `--config=testing`
to DUB. Furthermore, `fm-index` must be built:

```sh
cd external
make
make install
```

### Runtime Dependencies

The following additional software packages are required to run `dentist`
for testing:

- [SDSL][sdsl-lite]

Please see their own documentation for installation instructions.


[sdsl-lite]: https://github.com/simongog/sdsl-lite

Usage
-----

Use `snakemake --snakefile=Snakefile.testing`.
See `./snakemake/Snakefile.testing` for details.


Maintainer
----------

Arne Ludwig &lt;<arne.ludwig@posteo.de>&gt;


Contributing
------------

Contributions are warmly welcome. Just create an [issue][gh-issues] or [pull request][gh-pr] on GitHub. If you submit a pull request please make sure that:

- the code compiles on Linux using the current release of [dmd][dmd-download],
- your code is covered with unit tests (if feasible) and
- `dub test` runs successfully.


[gh-issues]: https://github.com/a-ludi/dentist/issues
[gh-pr]: https://github.com/a-ludi/dentist/pulls
[dmd-download]: https://dlang.org/download.html#dmd


License
-------

This project is licensed under MIT License (see license in [LICENSE](./LICENSE).

# caderno-zika

[![tests](https://github.com/gustavoM32/caderno-zika/actions/workflows/tests.yml/badge.svg)](https://github.com/gustavoM32/caderno-zika/actions/workflows/tests.yml)

Algorithms for competitive programming.

We're using [this tool](https://github.com/online-judge-tools/verification-helper) for testing. Supported judges are listed [here](https://online-judge-tools.github.io/verification-helper/document.html).

Ideally, for every algorithm `foo.cpp` there will be a corresponding file `foo.test.cpp` that uses it to solve some problem in an online judge. This helps finding any errors in the algorithm. The test file contains the link to the problem in the first line as `#define PROBLEM "link"`, followed by `#define IGNORE` if that problem is from an unsupported judge (see `example.test.cpp`).
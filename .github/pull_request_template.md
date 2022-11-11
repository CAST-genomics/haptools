Title your PR according to [the conventional commits specification](https://www.conventionalcommits.org/).

If your pull request resolves any open issues, please list them here at the top of the PR like this:

resolves #0

## Overview

Include a one sentence summary of the changes to our code and a separate paragraph explaining relevant motivation and context.

Also, list any dependencies that are required for this change.

## Usage and Documentation

You can call the new command like this...
Documentation for all of this is linked here...

## Implementation Details

The new command works by calling the following new classes...

### WickedAwesomeClass
The primary method for this class is...

### HelperForWickedClass
This works by...

## Tests

Please describe the tests that you ran to verify your changes. All tests should be added to our test suite in the `tests/` directory. Please also list any relevant details for your test configuration (ex: the version of pip you used, if it's relevant).

1. `test_a()`
    Handles the most basic case where...
2. `test_b()`
	Tests that we handle...

## Future work
I didn't have a chance to do...
In the future, we may want to refactor for...


## Checklist

* [ ] I have followed the [contributing guidelines](https://haptools.readthedocs.io/en/stable/project_info/contributing.html#how-to-fix-a-bug-or-implement-a-new-feature)
* [ ] I have adhered to the [style guidelines](https://haptools.readthedocs.io/en/stable/project_info/contributing.html#style)
* [ ] I have checked to ensure there aren't other open [pull requests](../../../pulls) for the same update/change
* [ ] (if applicable) I have added changes to the dependencies to the pyproject.toml file, and I've run `poetry lock` to ensure the lock file stays up to date

<!-- For an example that follows this PR outline, see [PR #43](https://github.com/CAST-genomics/haptools/pull/43) -->

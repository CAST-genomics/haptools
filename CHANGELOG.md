# Changelog

## [0.1.0](https://github.com/CAST-genomics/haptools/compare/v0.0.3...v0.1.0) (2023-01-31)


### Features

* a new `--no-normalize` parameter for `simphenotype` ([#156](https://github.com/CAST-genomics/haptools/issues/156)) ([24a0867](https://github.com/CAST-genomics/haptools/commit/24a0867ee8c4b1f495f5dde43600758daa2d916d))
* add `--seed` to `simphenotype` ([#162](https://github.com/CAST-genomics/haptools/issues/162)) ([9f7890a](https://github.com/CAST-genomics/haptools/commit/9f7890a6efd5be90717aaeeeab1f638512c6e2c4))
* change default verbosity to INFO ([#157](https://github.com/CAST-genomics/haptools/issues/157)) ([18ff090](https://github.com/CAST-genomics/haptools/commit/18ff090b431891b16b4aae9d297983d481a3d248))
* github actions for publishing ([#141](https://github.com/CAST-genomics/haptools/issues/141)) ([5cf6f29](https://github.com/CAST-genomics/haptools/commit/5cf6f29b64a49f9405ab78eebab5e0de20963376))
* Significantly decreased runtime of simgenotype  ([#163](https://github.com/CAST-genomics/haptools/issues/163)) ([de011d0](https://github.com/CAST-genomics/haptools/commit/de011d0106a4b50e0ea1e38eff6c567174a3a72f))
* updated `.vcf.gz` output in simgenotype ([#150](https://github.com/CAST-genomics/haptools/issues/150)) ([f61f613](https://github.com/CAST-genomics/haptools/commit/f61f613c55638efb610fb9ca03b0292f42dc2623))
* Updated simgenotype to allow only breakpoint generation ([#167](https://github.com/CAST-genomics/haptools/issues/167)) ([c0c3c97](https://github.com/CAST-genomics/haptools/commit/c0c3c97384c94056f79786f5297251bf9c29f2b9))


### Bug Fixes

* Added an additional test file containing only autosome centromeres ([#168](https://github.com/CAST-genomics/haptools/issues/168)) ([ca220e7](https://github.com/CAST-genomics/haptools/commit/ca220e7319f4b4d6594bad0042a9eb35c2a2d607))
* Added utilization of logging class to karyogram and simgenotype ([#159](https://github.com/CAST-genomics/haptools/issues/159)) ([bcf778a](https://github.com/CAST-genomics/haptools/commit/bcf778a1180d6208c7f24a6f63af1df830734dc5))
* case where genotypes are all the same but heritability is specified in `simphenotype` ([#145](https://github.com/CAST-genomics/haptools/issues/145)) ([063b411](https://github.com/CAST-genomics/haptools/commit/063b411302e4f299f02de82ad4040d2eb3fa8a61))
* Fixed logic error in simgenotype and docs errors ([#172](https://github.com/CAST-genomics/haptools/issues/172)) ([18dd27f](https://github.com/CAST-genomics/haptools/commit/18dd27fbe45c4af0081aa36928796cd61e30041e))
* issue where breakpoints weren't outputting ([#144](https://github.com/CAST-genomics/haptools/issues/144)) ([ff90bff](https://github.com/CAST-genomics/haptools/commit/ff90bff662b5b3ebc98278cf64d4919fbe4699c5))
* logging so that it doesn't affect the root logger ([#154](https://github.com/CAST-genomics/haptools/issues/154)) ([c1dc1ef](https://github.com/CAST-genomics/haptools/commit/c1dc1efed132f365c9268978be37bf25b4ef10b5))
* pin numpy to ensure a recent version of numpy ([#151](https://github.com/CAST-genomics/haptools/issues/151)) ([7c84a45](https://github.com/CAST-genomics/haptools/commit/7c84a454dfe3b57e559a3f4439b1f2dc390700d0))
* Update maps.rst ([#170](https://github.com/CAST-genomics/haptools/issues/170)) ([8ace893](https://github.com/CAST-genomics/haptools/commit/8ace8933644e98eda9483d47f4a895fe8cdd7831))


### Documentation

* better showcase example files and their locations in `simgenotype` ([#166](https://github.com/CAST-genomics/haptools/issues/166)) ([6995407](https://github.com/CAST-genomics/haptools/commit/69954075e547917882f72ca70776169254d16b99))
* updates for v0.1.0 ([#171](https://github.com/CAST-genomics/haptools/issues/171)) ([714142f](https://github.com/CAST-genomics/haptools/commit/714142f06204b284571f5b01c8a2bed3892dada7))

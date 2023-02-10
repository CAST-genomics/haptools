# Changelog

## [0.1.3](https://github.com/CAST-genomics/haptools/compare/v0.1.2...v0.1.3) (2023-02-10)


### Bug Fixes

* build workflow for comment bot ([#186](https://github.com/CAST-genomics/haptools/issues/186)) ([ebcf6cd](https://github.com/CAST-genomics/haptools/commit/ebcf6cd508701951ce2f1785eb3bcfee6aa69041))
* comment bot in build pipeline ([#184](https://github.com/CAST-genomics/haptools/issues/184)) ([b3fb839](https://github.com/CAST-genomics/haptools/commit/b3fb83939c6ec9141dc9c0d9b89836cfe0fee394))


### Documentation

* warn about pysam pip issue ([#187](https://github.com/CAST-genomics/haptools/issues/187)) ([bda0157](https://github.com/CAST-genomics/haptools/commit/bda015751dc5e5aa5d3672e77737d837cd116782))

## [0.1.2](https://github.com/CAST-genomics/haptools/compare/v0.1.1...v0.1.2) (2023-02-02)


### Bug Fixes

* add poetry readme to fix long_description when publishing to pypi ([#177](https://github.com/CAST-genomics/haptools/issues/177)) ([4050251](https://github.com/CAST-genomics/haptools/commit/405025161ef99b7c18c6adddb58d74bbcff93570))
* checkout in comment bot workflow for CD pipeline ([#181](https://github.com/CAST-genomics/haptools/issues/181)) ([9782d3d](https://github.com/CAST-genomics/haptools/commit/9782d3d5826a7029a74fb845dbfff29a7a0e8f2d))
* checkout ref in comment bot from CD pipeline ([#182](https://github.com/CAST-genomics/haptools/issues/182)) ([e3b92f6](https://github.com/CAST-genomics/haptools/commit/e3b92f6ba2764819d9936c7e37418bcd68469dfb))
* comment bot in CD pipeline ([#180](https://github.com/CAST-genomics/haptools/issues/180)) ([a2f66bd](https://github.com/CAST-genomics/haptools/commit/a2f66bd74dfb1fae8a81614912cba35366fe003b))
* indentation in checkout build ([#183](https://github.com/CAST-genomics/haptools/issues/183)) ([de431ab](https://github.com/CAST-genomics/haptools/commit/de431ab65a15f4dda2e7b6319608331b73ee5771))

## [0.1.1](https://github.com/CAST-genomics/haptools/compare/v0.1.0...v0.1.1) (2023-02-02)


### Documentation

* add dependency to conda env build ([#175](https://github.com/CAST-genomics/haptools/issues/175)) ([6626049](https://github.com/CAST-genomics/haptools/commit/66260498cbaf7fb652ababd6241cdfbc8a865dfd))

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

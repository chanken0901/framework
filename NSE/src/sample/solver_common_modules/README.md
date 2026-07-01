# Solver common modules

GPE と NSE の両方からそのまま `use` できる共通モジュール群です。

## Modules

- `mod_precision.f90`: `dp` の定義
- `mod_common_config.f90`: 格子・時間・出力・実行設定
- `mod_model_config.f90`: GPE/NSE 固有設定
- `mod_input_reader.f90`: `input.dat` の namelist 読み込み
- `mod_slf_output.f90`: 共通 field 出力 `*.slf` と `meta.json`

## Compile example

```bash
gfortran -c src/common/mod_precision.f90
gfortran -c src/common/mod_common_config.f90
gfortran -c src/common/mod_model_config.f90
gfortran -c src/common/mod_input_reader.f90
gfortran -c src/common/mod_slf_output.f90
```

## GPE usage

```fortran
use mod_common_config, only : simulation_config, init_simulation_config
use mod_model_config, only : gpe_config, init_gpe_config
use mod_input_reader, only : read_all_inputs
use mod_slf_output, only : write_field_slf
```

## NSE usage

```fortran
use mod_common_config, only : simulation_config, init_simulation_config
use mod_model_config, only : nse_config, init_nse_config
use mod_input_reader, only : read_all_inputs
use mod_slf_output, only : write_field_slf
```

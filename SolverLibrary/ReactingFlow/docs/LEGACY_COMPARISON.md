# 旧NSE反応流との直接比較

2026-10-01。式を別途書き直した模擬比較ではなく、旧NSEの実際のFortranソースを
検証専用ターゲットへ読み取り専用で組み込み、新ReactingFlowと比較する。
通常のReactingFlowライブラリ／実行ファイルにNSEへの依存は追加しない。

## 条件と結果

同じSI状態、等モル質量のA/B、NASA定比熱・生成エンタルピー、
Rusanov一次精度、SSPRK3、24セル、時間刻み0.1 µs、50ステップを使用。
非反応／発熱一段反応A→B、周期／両端鏡像の4組合せを比較した。
周期では組成正弦波、鏡像では圧力比2の衝撃波管を初期条件とする。
旧版の分子量kg/kmolと新版のkg/molは明示的に対応させる。

化学時間積分は旧版の陽的細分化と新版のDVODEで異なるため、ビット一致は要求しない。
両版のEuler時間積分と旧化学積分の実コードを使用して、同じStrang順序で比較する。

| 条件 | 全保存変数の最大規格化差 |
|---|---:|
| 非反応・周期 | 4.35e-15 |
| 非反応・鏡像 | 1.06e-10 |
| 反応・周期 | 5.39e-10 |
| 反応・鏡像 | 1.08e-9 |

規格化差はセル・変数ごとの `abs(new-old)/max(1,abs(new))` の最大値。
全ケースで1e-7未満、質量と境界寄与込みのエネルギー保存も合格。
この基準はSI成分ごとの数値比較であり、任意単位系共通の誤差ノルムではない。

## 実行

FrameWorkルートで、既存のFortranビルド環境（CMake、Ninja、gfortran）を使用する。
DVODEは通常ビルドと同じ固定リビジョン。オフライン時は既存手順に従って
`RF_DVODE_SOURCE_DIR` を追加指定する。

PowerShell:

```powershell
cmake -S SolverLibrary/ReactingFlow -B build/reactingflow-legacy-check -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_COMPILER=gfortran -DRF_LEGACY_NSE_SOURCE_DIR="$PWD/SolverLibrary/NSE"
cmake --build build/reactingflow-legacy-check --target rf_legacy_unit
ctest --test-dir build/reactingflow-legacy-check -R rf_legacy_unit -V
```

Linux:

```bash
cmake -S SolverLibrary/ReactingFlow -B build/reactingflow-legacy-check -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_COMPILER=gfortran -DRF_LEGACY_NSE_SOURCE_DIR="$PWD/SolverLibrary/NSE"
cmake --build build/reactingflow-legacy-check --target rf_legacy_unit
ctest --test-dir build/reactingflow-legacy-check -R rf_legacy_unit -V
```

`RF_LEGACY_NSE_SOURCE_DIR` を指定しない通常ビルドでは追加ターゲットを作らない。
旧ソースが存在しない場合は比較を偽装せずCMake構成時に失敗する。

## 限界

1D断面相当の非粘性・一段反応の対応比較であり、旧CFDの全機能同値ではない。
旧MUSCLとの比較ではなく、両者が持つ一次精度Rusanovを対応させた。
分子輸送、詳細機構、一般座標・ノズル、多次元・並列・GPUの移行比較は別途必要。

# GPE SolverLibrary

`gp3d`は三次元Gross-Pitaevskii方程式をsplit-operator擬スペクトル法で解く
モジュール型ソルバーです。

- CPU逐次: 自前DFT / FFTW
- CPU MPI: 従来のzスラブ分割、または2DECOMP&FFTによるペンシル分割FFT
- 単一GPU: CUDA / cuFFT
- 複数GPU: MPI / cuFFTMp（スラブ分割またはペンシル分割）

実行環境へコピーするファイルは`gp3d/solver_manifest.yaml`のプロファイルで決まります。
ケース依存条件はSolverLibraryへ書かず、研究プロジェクトの`cases/<case_id>/case.yaml`へ
記録してください。

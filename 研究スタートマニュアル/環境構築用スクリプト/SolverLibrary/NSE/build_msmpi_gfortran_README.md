# MS-MPI + gfortran ビルドスクリプト

## 対象ファイル

このスクリプトは、以下の2つのFortranファイルを同じフォルダに置いて使います。

```text
module_mpi_TW_NEC.f90
3d_mpi_openmp_ver2_TGV_KEEP.f90
build_msmpi_gfortran.bat
```

## 実行方法

cmdで対象フォルダに移動して、

```bat
build_msmpi_gfortran.bat
```

を実行します。

成功すると、

```text
solver.exe
```

が生成されます。

## MPI実行

例：

```bat
mpiexec -n 4 solver.exe
```

## 入れている重要オプション

```text
-fallow-invalid-boz
```

MS-MPI の `mpif.h` で出る BOZ literal のエラー対策です。

```text
-fallow-argument-mismatch
```

`include 'mpif.h'` を用いた古いMPIコードで、REAL(4), REAL(8), INTEGER などの引数型チェックが厳しくなって止まる問題への対策です。

```text
-lmsmpifec -lmsmpi
```

MS-MPI の Fortran MPI 呼び出しをリンクするために必要です。
`-lmsmpi` だけでは `mpi_isend_` などが未解決になる場合があります。

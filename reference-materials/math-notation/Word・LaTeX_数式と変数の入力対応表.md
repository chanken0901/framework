# Word・LaTeX 数式と変数の入力対応表

作成日：2026-11-10

**Excel版**：[Word・LaTeX 数式と変数の対応表（.xlsx）](../../outputs/math-notation-01a1bdee/Word_LaTeX_数式と変数の対応表.xlsx)
を利用してください。下記の内容に加え、流体の保存則・物性式など指定の9式、変数の意味・単位・適用条件を収録しています。このMarkdownは変換前の参考資料として残しています。

特定のフレームワークに依存しない、論文・レポート用の数式入力資料です。
「表示例」は完成した数式、「入力例」はキーボードで入力する文字です。
Markdownの数式表示に対応したプレビューで読むと、表示例が数式として見えます。

## 1. 入力を始める前に

### Word

通常の本文ではなく、**数式欄**に入力します。

1. 「挿入」→「数式」で数式欄を作ります。Windowsでは`Alt`＋`=`でも作成できます。
2. 数式用タブの入力方式で「UnicodeMath（Unicode）」または「LaTeX」を選びます。
3. 以下の表の、その方式に対応する文字列を入力します。
4. 「変換」→「2次元形式／Professional」で、分数などの見た目へ変換します。表示名は版・言語によって異なります。

UnicodeMathでは、`\alpha`などの記号コマンドの直後に半角Spaceを押すと変換されます。
アクセントなどはSpaceを2回必要とする場合があります。以下の`〈Space〉`は説明用であり、その文字自体は入力しません。
WordのLaTeXモードでは、表のLaTeX列の**数式部分だけ**を数式欄へ入力します。
Wordのバージョンによって対応範囲が異なるため、複雑な式は変換後の見た目も照合してください。
操作の根拠：[Microsoftの数式入力説明](https://support.microsoft.com/en-us/word/linear-format-equations-using-unicodemath-and-latex-in-word)。

### LaTeX

LaTeX文書では数式を区切り記号で囲みます。表には、この区切りの**内側に書く部分**を示しています。

| 使い方 | 入力例 |
|---|---|
| 文章中の数式 | `速度は \(u_i\) とする。` |
| 独立した行の数式 | `\[ E=mc^2 \]` |
| 番号付きの数式 | `\begin{equation} E=mc^2 \end{equation}` |

行列・複数行の式などを使う文書では、本文開始前に`\usepackage{amsmath,amssymb}`を記述します。
Wordの数式欄へは`\documentclass`や`\usepackage`を書きません。
参考：[Overleafの数式モード説明](https://www.overleaf.com/learn/latex/Mathematical_expressions)。

## 2. 添字・分数・基本演算

Word欄は**UnicodeMathモード**の例です。WordでLaTeXモードを選んだ場合は右端の列を使います。

| 意味 | 表示例 | Word：UnicodeMath | LaTeX／Word：LaTeX |
|---|---|---|---|
| 変数 | $x$ | `x` | `x` |
| 下付き添字 | $u_i$ | `u_i` | `u_i` |
| 複数文字の添字 | $u_{ij}$ | `u_(ij)` | `u_{ij}` |
| 上付き・べき乗 | $x^2$ | `x^2` | `x^2` |
| 複数文字の上付き | $x^{n+1}$ | `x^(n+1)` | `x^{n+1}` |
| 上下の添字 | $u_i^{n+1}$ | `u_i^(n+1)` | `u_i^{n+1}` |
| 分数 | $\frac{a}{b}$ | `a/b` | `\frac{a}{b}` |
| 式を含む分数 | $\frac{a+b}{c+d}$ | `(a+b)/(c+d)` | `\frac{a+b}{c+d}` |
| 平方根 | $\sqrt{x}$ | `\sqrt(x)` | `\sqrt{x}` |
| 大きさが変わる括弧 | $\left(\frac{a}{b}\right)$ | `(a/b)` | `\left(\frac{a}{b}\right)` |
| 掛け算の記号 | $a\times b$ | `a \times b` | `a\times b` |
| 中点 | $a\cdot b$ | `a \cdot b` | `a\cdot b` |
| プラスマイナス | $a\pm b$ | `a \pm b` | `a\pm b` |
| 等しくない | $a\ne b$ | `a \ne b` | `a\ne b` |
| 以下・以上 | $a\le b,\ a\ge c$ | `a \le b, a \ge c` | `a\le b,\ a\ge c` |
| ほぼ等しい | $a\approx b$ | `a \approx b` | `a\approx b` |
| 比例する | $a\propto b$ | `a \propto b` | `a\propto b` |
| 無限大 | $\infty$ | `\infty` | `\infty` |

LaTeXの`{...}`は「ここまでをひとまとまりにする」ための記号で、通常は表示されません。
例えば`x^n+1`は$x^n+1$、`x^{n+1}`は$x^{n+1}$で、意味が異なります。
WordのUnicodeMathでは添字のまとまりに丸括弧を使う点が大きな違いです。

## 3. ギリシャ文字とよく使う記号

この表の入力はWordの両方式とLaTeXで共通です。UnicodeMathでは記号名の直後にSpaceを押して変換します。
「用途例」は一例です。記号の意味は論文ごとに定義してください。

| 表示 | 読み方 | 入力 | 用途例 |
|---|---|---|---|
| $\alpha$ | アルファ | `\alpha` | 係数 |
| $\beta$ | ベータ | `\beta` | 係数 |
| $\gamma$ | ガンマ | `\gamma` | 比熱比など |
| $\delta$ | デルタ | `\delta` | 微小変化、デルタ記号 |
| $\Delta$ | 大文字デルタ | `\Delta` | 差、変化量 |
| $\epsilon$ | イプシロン | `\epsilon` | 小さい量、誤差 |
| $\varepsilon$ | イプシロンの異体字 | `\varepsilon` | 散逸率など |
| $\theta$ | シータ | `\theta` | 角度 |
| $\kappa$ | カッパ | `\kappa` | 係数 |
| $\lambda$ | ラムダ | `\lambda` | 波長、固有値 |
| $\mu$ | ミュー | `\mu` | 粘性係数など |
| $\nu$ | ニュー | `\nu` | 動粘性係数など |
| $\pi$ | パイ | `\pi` | 円周率 |
| $\rho$ | ロー | `\rho` | 密度 |
| $\sigma$ | シグマ | `\sigma` | 標準偏差、応力など |
| $\tau$ | タウ | `\tau` | 時間尺度、応力など |
| $\phi$ | ファイ | `\phi` | 位相、ポテンシャル |
| $\varphi$ | ファイの異体字 | `\varphi` | 位相など |
| $\psi$ | プサイ | `\psi` | 波動関数など |
| $\omega$ | オメガ | `\omega` | 角周波数など |
| $\Omega$ | 大文字オメガ | `\Omega` | 領域、角速度など |
| $\partial$ | 偏微分記号 | `\partial` | 偏微分 |
| $\nabla$ | ナブラ | `\nabla` | 勾配、発散、回転 |
| $\hbar$ | エイチバー | `\hbar` | 換算プランク定数 |

英字の`v`とギリシャ文字の$\nu$、英字の`p`と$\rho$などを混同しないようにします。
大文字の`\Sigma`は文字$\Sigma$、`\sum`は総和演算子$\sum$です。
参考：[Overleafの記号一覧](https://www.overleaf.com/learn/latex/List_of_Greek_letters_and_math_symbols)。

## 4. 微分・積分・和

WordのUnicodeMath例には、貼り付けやすいように`∂`・`∫`・`∑`などの記号を直接使っています。
それぞれ`\partial`・`\int`・`\sum`＋Spaceで入力してから式を続ける方法でも構いません。

| 意味 | 表示例 | Word：UnicodeMath | LaTeX／Word：LaTeX |
|---|---|---|---|
| 常微分 | $\frac{df}{dx}$ | `df/dx` | `\frac{df}{dx}` |
| 2階微分 | $\frac{d^2f}{dx^2}$ | `(d^2 f)/(dx^2)` | `\frac{d^2f}{dx^2}` |
| 偏微分 | $\frac{\partial f}{\partial x}$ | `(∂f)/(∂x)` | `\frac{\partial f}{\partial x}` |
| 2階偏微分 | $\frac{\partial^2 f}{\partial x^2}$ | `(∂^2 f)/(∂x^2)` | `\frac{\partial^2 f}{\partial x^2}` |
| 混合偏微分 | $\frac{\partial^2 f}{\partial x\partial y}$ | `(∂^2 f)/(∂x ∂y)` | `\frac{\partial^2 f}{\partial x\partial y}` |
| ラプラシアン | $\nabla^2 f$ | `∇^2 f` | `\nabla^2 f` |
| 不定積分 | $\int f(x)\,dx$ | `∫ f(x) dx` | `\int f(x)\,dx` |
| 定積分 | $\int_a^b f(x)\,dx$ | `∫_a^b f(x) dx` | `\int_a^b f(x)\,dx` |
| 閉曲線・閉曲面の積分記号 | $\oint_C f\,ds$ | `∮_C f ds` | `\oint_C f\,ds` |
| 総和 | $\sum_{i=1}^{N} a_i$ | `∑_(i=1)^N a_i` | `\sum_{i=1}^{N} a_i` |
| 総積 | $\prod_{i=1}^{N} a_i$ | `∏_(i=1)^N a_i` | `\prod_{i=1}^{N} a_i` |

積分の`\,`は小さな空白です。掛け算の意味を追加する記号ではありません。
微分の$d$を立体で組む流儀では、LaTeXの`d`を`\mathrm{d}`に置き換えます。
例：`\frac{\mathrm{d}f}{\mathrm{d}x}`、`\int_a^b f(x)\,\mathrm{d}x`。
投稿先の指定に合わせ、同じ文書で表記を統一します。

## 5. ベクトル・平均・文字の書体

ここは**LaTeX、またはWordのLaTeXモード**を使うと対応をそろえやすい部分です。
WordのUnicodeMathでアクセントを付ける場合、例えば矢印付きの$u$は`u\vec`に続けてSpaceを押して変換します。

| 意味・区別 | 表示例 | LaTeX／Word：LaTeX |
|---|---|---|
| スカラー変数 | $u$ | `u` |
| 矢印付きベクトル | $\vec{u}$ | `\vec{u}` |
| 太字ベクトル | $\mathbf{u}$ | `\mathbf{u}` |
| ベクトルの成分 | $u_i$ | `u_i` |
| 行列を太字で表す例 | $\mathbf{A}$ | `\mathbf{A}` |
| バー付き文字 | $\bar{u}$ | `\bar{u}` |
| ハット付き文字 | $\hat{u}$ | `\hat{u}` |
| チルダ付き文字 | $\tilde{u}$ | `\tilde{u}` |
| ドット付き文字 | $\dot{x}$ | `\dot{x}` |
| 二つのドット | $\ddot{x}$ | `\ddot{x}` |
| 平均の山括弧 | $\langle u\rangle$ | `\langle u\rangle` |
| 内積 | $\mathbf{a}\cdot\mathbf{b}$ | `\mathbf{a}\cdot\mathbf{b}` |
| 外積（クロス積） | $\mathbf{a}\times\mathbf{b}$ | `\mathbf{a}\times\mathbf{b}` |
| 勾配 | $\nabla f$ | `\nabla f` |
| 発散 | $\nabla\cdot\mathbf{u}$ | `\nabla\cdot\mathbf{u}` |
| 回転 | $\nabla\times\mathbf{u}$ | `\nabla\times\mathbf{u}` |
| 絶対値 | $\lvert x\rvert$ | `\lvert x\rvert` |
| ノルム | $\lVert\mathbf{u}\rVert$ | `\lVert\mathbf{u}\rVert` |
| 説明用の添字 | $u_{\mathrm{ref}}$ | `u_{\mathrm{ref}}` |
| 単位を立体にする | $10\,\mathrm{m}/\mathrm{s}$ | `10\,\mathrm{m}/\mathrm{s}` |
| 標準の関数名 | $\sin x,\ \exp x,\ \log x$ | `\sin x,\ \exp x,\ \log x` |

変数は通常斜体、単位・説明用の添字・標準関数名は立体にします。
例えば`sin x`と入力すると$s i n x$のように変数の積と区別しにくくなります。関数には`\sin x`を使います。
`\mathbf`は主に英数字の太字向けです。ギリシャ文字を太字にするLaTeXでは、`amsmath`の`\boldsymbol{\omega}`などを使います。
Word側の対応や書体は版に依存するため、変換結果を確認してください。
参考：[Microsoft 365のLaTeX対応一覧](https://learn.microsoft.com/en-us/office/math/latex)、[Overleafの数式書体](https://www.overleaf.com/learn/latex/Mathematical_fonts)。

バー・ハット・チルダは記号だけでは意味が決まりません。
平均、正規化、フーリエ変換など、何を表すかを本文で定義します。

## 6. 方程式全体の入力例

### 二次方程式の解

$$x=\frac{-b\pm\sqrt{b^2-4ac}}{2a}$$

Word・UnicodeMath：

```text
x=(-b±√(b^2-4ac))/(2a)
```

LaTeX／Word・LaTeX：

```latex
x=\frac{-b\pm\sqrt{b^2-4ac}}{2a}
```

### 1次元の拡散方程式

$$\frac{\partial u}{\partial t}=D\frac{\partial^2u}{\partial x^2}$$

Word・UnicodeMath：

```text
(∂u)/(∂t)=D (∂^2 u)/(∂x^2)
```

LaTeX／Word・LaTeX：

```latex
\frac{\partial u}{\partial t}=D\frac{\partial^2u}{\partial x^2}
```

### 連続の式の表記例

$$\frac{\partial\rho}{\partial t}+\nabla\cdot(\rho\mathbf{u})=0$$

LaTeX／Word・LaTeX：

```latex
\frac{\partial\rho}{\partial t}+\nabla\cdot(\rho\mathbf{u})=0
```

この例では$\rho$を密度、$\mathbf{u}$を速度ベクトル、$t$を時刻と定義します。
ベクトル全体の$\mathbf{u}$と、成分を表す$u_i$は書き分けます。

## 7. 行列・複数行の式

### 2行2列の行列

$$A=\begin{pmatrix}a&b\\c&d\end{pmatrix}$$

Word・UnicodeMathでは、次を入力して変換します。

```text
A=(\matrix(a&b@c&d))
```

`&`は列の区切り、`@`は行の区切りです。
LaTeXでは次を使います（`amsmath`を読み込んでおきます）。

```latex
A=\begin{pmatrix}a&b\\c&d\end{pmatrix}
```

WordのLaTeX取り込みで行列環境が使えるかは版によって異なります。
変換されない場合は、Wordの「行列」「かっこ」のテンプレートを使うか、UnicodeMathの例へ切り替えます。

### 等号の位置をそろえる

LaTeX（`amsmath`を使用）：

```latex
\begin{align}
a &= b+c \\
  &= d+e
\end{align}
```

`&`が位置をそろえる基準、`\\`が改行です。この`align`環境の外側をさらに`\[...\]`で囲みません。
Word・UnicodeMathでは、数式欄に次を入力して変換できます。

```text
\eqarray(a&=b+c@&=d+e)
```

WordはLaTeX文書全体を組版するソフトではありません。LaTeX側の環境・パッケージ・独自マクロがそのまま全部使えるわけではありません。
参考：[MicrosoftのUnicodeMath例](https://support.microsoft.com/en-us/word/linear-format-equations-using-unicodemath-and-latex-in-word)、
[Overleafの式の整列](https://www.overleaf.com/learn/latex/Aligning_equations_with_amsmath)。

## 8. 入力時の確認事項

- コマンドは半角で入力します。`\`はバックスラッシュです。日本語環境では同じ文字が円記号のように表示される場合があります。
- 表のコードを囲むバッククォートは入力しません。
- Wordで変換されない場合は「数式欄か」「入力方式が合っているか」を最初に確認します。
- LaTeXは大文字と小文字を区別します。`\delta`と`\Delta`は別の記号です。
- プログラムの掛け算`*`は、論文中では必ずしも`*`で表示しません。通常の積は$ab$、記号が必要なら$\times$や$\cdot$などを使い分けます。
- 方程式の入力後は、添字の範囲、分母の範囲、括弧、ベクトルの太字、微分の階数を確認します。
- 本表は公式構文資料を参照した入力例です。個々のWord環境で全例を実操作検証したものではありません。

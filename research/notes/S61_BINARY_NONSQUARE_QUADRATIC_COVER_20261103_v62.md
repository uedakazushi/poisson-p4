# K=0 三 core 開で六二次式は全点で独立

2026-11-03 / v62。全 S61 の元135式から既に構成された、七 light 変数における六つの正確な重み二次残余を用いる。結論は

    U={aFBC≠0} の全複素点で rank(6×28 二次係数行列)=6

である。ここで

    F=500a+90e²+3e−24,
    B=(50a+18e²−9e−2)d+60e−40,
    C=100a+18e²−9e−2.

従来は階数欠損集合が曲線以下という結論だった。今回はその集合自体が U 上で空であることを QQ 上の exact ideal 計算で示した。重み三・四の欠損集合を空とする結果ではない。

## 元式と可逆な三 pivot

入力は `reports/S61_BINARY_HYPERPLANE_QUADRATIC_CODIMENSION_20261102.json` の六行であり、その全入力ハッシュを照合する。この旧 producer は元135式から陰関数の二次補正を再構成し、各行の content が a,F,B と有理定数だけからなることを確認している。

七 light 変数の二次単項式28個を列に並べ、旧報告の三 pivot を正確に再実行する。各 pivot の分子・分母が a,F,B,C の単元であることを多項式完全除算で再確認する。残る三行の既知単元分母を払い、25列の残余行列 N を得る。行の追加 content は有理定数だけであり、未知 core 因子を捨てていない。

したがって U 上では元の階数は `3+rank N` である。

## 全 residual minor ideal

N の全3×3小行列式を計算する。恒等零を除き、a,F,B,C の因子だけを完全に除いて有理定数で monic に正規化する。同じ多項式を重複排除すると1010個になる。この除去は U 上で可逆である。

新しい変数 z を取り、1010多項式と `z*a*F*B*C−1` が生成する QQ[a,e,d,z] の ideal J を作る。Macaulay2 1.24.11 の正確な Gröbner 計算は

    gens gb J = matrix {{1}},    dim(R/J)=−1

を返す。producer は結果の文字列だけでなく `assert(1%gb J==0)` も実行する。従って J は単位 ideal であり、U 上の N の階数欠損点はない。有限体の非零標本や gcd=1 から単位 ideal を推測した証明ではない。

この exact ideal 計算は任意の標数零係数体への拡大後にも成立するので、全複素 core 点を扱う。

## 再現と残る課題

    source /workspace/cloud-setup/activate.sh
    .venv/bin/python code/s61_binary_nonsquare_quadratic_cover_20261103_v62.py \
      --output /tmp/nonsquare-quadratic-cover-new

実行時間21.31秒。凍結した [report](../../reports/S61_BINARY_NONSQUARE_QUADRATIC_COVER_20261103_v62/report.json) に全 minor の選択列、除いた既知因子の指数、入力ハッシュ、生成 M2 入力・出力のハッシュを保存した。M2 入力は同ディレクトリの `cover.m2`、結果は `cover.log`。

全 marked 非平方点の除去には、[残余 core の零次元化](S61_NONSQUARE_MARKED_REMAINING_20261103.md) が必要である。二次階数については、U 上の未処理曲線・点はなくなった。今後は二次式の一次倍、heavy 三次ブロック、全三次、および四次全射の欠損を扱う。この補題単独から横断長65や周囲の t 冪零性を主張しない。

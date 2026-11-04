# NEW25 の七つの低次数不変式：全係数展開の独立検証

2026-11-04。既受理の実際の source-five 回路1032 nodes を入力とする。
著者の native worker や producer は再実行しない。beta16 の祖先を別に有理
算術で評価し、c=1、u=aX の下で `1080BD*beta_i` の全16多項式を作る。
その後は **Python 任意精度整数・疎係数辞書**による独立の word 演算を使う。

保存された native input の全16 source 多項式、136 monomial の順序、全32
covariant 係数行も全て照合する。指数を base256 に詰めるが、全中間 word の
beta 次数は12以下、元の各変数の次数は4以下なので carry は起こらない。

F4 は moment square、p0,…,p4 は determinant pencil、N12 は定義の trace
word から新たに作る。pencil には著者の permutation determinant と違う
traceless4×4 Newton identity を使う。N12 の六次 word は逆の cyclic 順序
`tr(D²CD²C)` で計算する。native の整数出力だけを正しいと仮定しない。

全七多項式の各係数を、整数と B の冪を払って照合する。D の冪は両辺とも
4、8、12のまま保持する。X の parity から `Y=X²` への写像、全係数の整数
content、B の最大 monomial cancellation、c の復元冪2、4、6も検査する。
term 数は F4=96、p0=336、p1=1040、p2=1061、p3=973、p4=893、N12=1700。
従って source-four形状 `(a=u/X,Y=X²,W,B)` 上の七つの完全な係数 map を
この範囲で受理する。新しい beta 点は無く、c=1 は記号的な斉次化の操作。

## 後続の degree24 relation に必要な別の gate

F4,p1,…,p4,N12 の weighted degree24 monomial は41個である。この中で
有限体の kernel が見つかっても、それだけでは方程式にならない。非零な
QQ係数候補について**全 source 係数の厳密な零**を別に確認する必要がある。
さらに全 W16 上で選んだ六不変式の代数的独立性があれば、その候補は普遍的
恒等式ではなく非零な beta 多項式となる。source-five が実際の NEW quotient
を支配する gate と、既受理の低次数除外 gate を合わせて初めて NEW の定義式
次数を結論できる。degree24 の全 Inv24 空間と、この41部分空間は区別する。

本検証はその41 relation test を実行せず、NEW 方程式や全 Inv24 の除外を
主張しない。有限指数が未計算の source quotient と実際の quotient も
同一視しない。支配性による恒等式の下降には、その指数1は不要である。

```sh
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/q_new_family_low_invariant_review_20261104_v63.py --output work/new-family-low-invariant-independent
```

60秒・2GB上限。係数辞書の term cap は150000。高次数への自動拡張は無い。

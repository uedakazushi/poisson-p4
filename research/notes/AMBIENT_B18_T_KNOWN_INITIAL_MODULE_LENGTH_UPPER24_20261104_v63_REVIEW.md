# 独立監査：実 T の長さ区間 10–24

**受理する結論は `10 <= length(T) <= 24` である。** 既知の homogeneous
initial 部分加群が真の initial 加群全体に等しいとは主張しない。元の
producer は実行せず、source/code/note/evidence hash を固定し、別の有理数
多項式 engine で検算する。

## 多項式と source の監査

実五生成825関係の frame は、八個の actual CA-F9 parameter
`(zA0,...,zA6,rA)` による R0/(ell) 上のものを保つ。元37一次初期関係を
保存された source から再構成する。二つの追加 quadratic は、元178本の
うちの実 source quadratic へ戻る直接多項式恒等式で確認する。

```
e4*z5²,
e3*r²+(11767/3652)*e4*z5*r.
```

178本中の177本は、実一次関係の変数倍および full Schur の literal d3
二次関係の真の線形結合である。最後の canonical rank1 representative
は、実620 source quadratic との差がその177空間に属するため、同じく
真の quadratic initial 関係である。この source bridge と actual HF(5,3,2)
は独立監査済み packet の hash に結ぶ。

本監査の polynomial key は八指数と最後の target成分である。元producerの
最初の適用可能除数を取る方針と反対に、最後の適用可能除数を取って厳密な
Fraction arithmetic で除算する。元の178本全てがこの39本へ零簡約される
ことと、同じ target の全 S-pair が零簡約されることを独立に確認する。
正確な数は **39 generators と133同一target S-pairs** である。

従ってこの39本は、指定 position-over-lex module 順序の全次数 Gröbner
基底である。先頭項から、target0,1,2では定数だけ、target3では定数と r、
target4では `r^d` と `z5*r^(d-1)` が残る。既知部分加群の商の Hilbert 級数は

```
5+3*t+2*t²/(1-t).
```

これは全次数の Buchberger 判定であり、有限次数で見た pattern の外挿
ではない。

## 実 annihilator と filtration

既知の actual T は、正しい比較完全列から IA*IB の **二乗**で零化される。
独立監査済み348基底は同じ ell=0 actual conductor pair を全次数の単位
変換で正規化し、その商の最大 ideal 十乗零を証明している。単位による
実局所座標変換は最大 ideal を保つので、これはこの五生成表示の Rbar-adic
filtration に対する `m^10*T=0` である。単なる type 同定や、H20上の naïve
kernel の結果としては用いない。

既知 homogeneous 部分加群の商から gr_m(T) への全射がある。従って各層の
次元は `5,3,2,2,...` 以下であり、degree0から9だけを足せば

```
length(T) <= 5+3+8*2 = 24.
```

独立確認済みの最初三層から下界10がある。既証の
`mu_D K=length(T)-4` と `rank_D K=4` により、最小生成数は6–20、DVR torsion
部分の最小生成数は2–16となる。後者は torsion exponent、torsion長、または
H0(D0) 長の上界ではない。新たな HF3、正確な T長、nil 三乗の結論はない。

今回の監査は非斉次 annihilator の式をそのまま gr(T) に代入しない。
さらに小さい上界のためには、別の source-bound initial annihilator gate
が必要であり、この24の受理範囲へ混ぜない。

## 再現

```
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/ambient_B18_T_known_initial_module_length_upper24_review_20261104_v63.py
```

新しい source coefficient、jet、rank search、348 scalar matrix は不要。
60秒・2GB・32MiB・150000 trace/係数項の固定 guard 内で検証する。

# L-frame からの実319 source の新しい方向

既存の structured reduction は canonical11 の generic rank10 を証明した。
その凍結済み入力を読むと、genuine relative 行 G6 と G1 は残り31係数で
既に同一である。full40では差が `L²*(L,x0,y1,y2,y3)` の五係数に限られる。
両行の同じ clearing は `delta=-64000*C²`。

正確な元八変数・全91係数で次を検証する。

```
delta*(G6-G1)-8100*q1*E = Z*A+W*B,
Z=3*v0-v1+v6,  W=v3.
```

E は linear、A,B は重み2の29 monomial（H を含む）の多項式である。
W を先に取り除き、v6 の多項式を v1-3v0 で割る literal monomial formula
で A,B を得る。solver・新 pivot・新点は用いない。

六つの実二次 source の正確な変換は

```
a*ell*Z=(10a+3e²-3e)q1-3(e-1)q3-10a*q4,
ell*W=-(de+8)q1+8q2+dq3.
```

従って前の恒等式に a*ell を掛けると、二次四行と cubic 二行だけを使う
実 full319 source-kernel vector を得る。genuine normalized43,49
（元135の実 cubic20,61）の cofactor はそれぞれ `-a*delta*ell`,
`a*delta*ell`。全33 quartic cofactor は零。
全239の重み4係数を直接検証し、実319順の vector を保存する。
全ての新しい cofactor は QQ[a,e,d] の多項式であり、未知因子は反転しない。

この方向は FULLU 上で旧37方向から独立である。既存 direct Bianchi15 の
quartic projection は既知単元11 pivotを持つ。quartic-zero の組合せを
その residual4へ移す。先の三つの cubic pivot は G18*v2, G18*v3,
G86*v0 であり、新 vector は全て零なのでその三係数は零になる。
最後の residual は G20*v1 と G20*v6 の係数和が零である。一方、新 vector は

```
coefficient(G20*v1)+coefficient(G20*v6)
  = -1024000*a*C³*J/243
```

という FULLU の既知単元を持つ。従って新 vector の係数も零。
旧 quadratic22 は cubic・quartic係数を持たないためこの議論に寄与しない。
これで全点独立方向は38、80まで残り42である。これは新しい一次方向一つの
証明であり、full239・全kernel尽くし・全cubic57を証明しない。
独立レビュー完了までは primary 証明として扱う。

```sh
source /workspace/cloud-setup/activate.sh
.venv/bin/python code/s61_nonsquare_ell_compact_relation_20261104_v63.py
```

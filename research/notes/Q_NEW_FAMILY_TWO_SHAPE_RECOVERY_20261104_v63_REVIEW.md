# q/B 回収と二 shape 制約の独立検証

2026-11-04。対象は `Q_NEW_FAMILY_MOBIUS_SHAPE_RECOVERY_20261104_v63` と
`Q_NEW_FAMILY_TWO_SHAPE_HOMOGENIZATION_20261104_v63`。
受理済み七写像の全係数を入力にし、別の sparse Fraction 多項式演算で検証する。
原著 producer、失敗した generic rational cancellation は再実行しない。

## 有理逆式と分母

既知 involution の固定 source 体は `QQ(xi,r,s,q,B)` であり、
`h=s²−25rq²`、`F=100r(r−1)q²−4rs+s²+4r`、`alpha=2/405` と置く。
物理不変式は `F4=alpha*xi²*F`、
`p0=alpha²*xi⁴*B*D/h`、
`D=2359296000r⁴(r−1)q²(5q−s+1)` である。

`Z0=p0/F4²`、`R=(p3+4p4)/p0` とすると、全係数恒等式として

```text
R=(80qr²−6r²s+6r²−11rs+10r+s)/(8r(r−1)(5q−s+1)),
T=(r−1)R−2r,
L=8Rr(r−1)−6r²−11r+1,
Nq=(s−1)L−(r−1),
q=Nq/(40rT),
B=Z0*h*F²/D.
```

Möbius determinant は `−40r(r−1)²*((10r−1)s−10r)` で非零。
両向き合成を分子の係数恒等式で確認し、従って**この source 体について**
`QQ(xi,r,s,q,B)=QQ(xi,r,s,R,Z0)` が成立する。
実 NEW target 体の index1 の主張ではない。

`J=(10r−1)s−10r` とし、

```text
Hn=64rT²s²−Nq²,
Fn=(r−1)Nq²+16rT²(s²−4rs+4r).
```

すると `h(q)=Hn/(64rT²)`、`F(q)=Fn/(16rT²)`、
`5q−s+1=(r−1)J/(8rT)`、
`D(q)=184320r(r−1)²Nq²J/T³`。
表示する dense chart は
`xi*r*s*Z0*(r−1)*T*Nq*J*Hn*Fn != 0` であり、各 factor の非零は
記号多項式で確認する。標本点による仮定ではない。

## 六有理係数と三制約

七つの受理済み map の `xi` 因子を除いた分子を
`p_i=(A_i+B*L_i)/h`、`N12=(B*M0+B²*M1)/h²` と書く。
この節では `F4=alpha F` とし、必要な `xi` 因子は正規化で相殺される。
直接 B を代入すれば

```text
Z_i=p_i/F4²=A_i/(alpha²*h*F²)+Z0*L_i/(alpha²*D), i=1,2,
ZN=N12/F4³=Z0*M0/(alpha³*h*F*D)+Z0²*F*M1/(alpha³*D²).
```

別の sparse 演算で `q=Nq/(40rT)` を全六分子へ homogenize し、
保存された530,105,554,109,807,205項を全係数比較する。
六つの reduced rational pair は、gcd を仮定せず named factor の交差積で検証する。
三つの cleared constraint も全係数を再構築した。

| 制約 | shape `(r,s)` bidegree 上界 | R degree 上界 | 係数項数合計 |
|---|---:|---:|---:|
| Z1 | (9,3) | 3 | 284 |
| Z2 | (14,5) | 5 | 1076 |
| ZN | (10,4) | 4 | 471 |

これらは表示した chart 上の厳密な同値制約である。分母零上では余分な成分が
生じ得るので、global extension・properness・elimination・最小次数・実定義方程式を
主張しない。未知の shape は `(r,s)` の二つに減るが、その消去や各解の次数、
実 target 体との有限 index は未解決である。

元の generic SymPy cancellation が60秒で停止した記録は保存する。
新点、source26、resultant、GB、generic multivariate gcd は用いない。

再現:

```sh
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/q_new_family_two_shape_recovery_review_20261104_v63.py --output work/v62-review/two-shape-recovery-independent
```

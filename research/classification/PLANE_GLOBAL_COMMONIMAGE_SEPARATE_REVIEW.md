# 全重複度 common-image 四全源の独立監査

2026-11-07。主義務 G。別作者による幾何証明の読解と、既存の独立 Fraction
ライブラリだけを用いる全行再構成。producer/CAS は import しない。
主原稿・global 台帳・Git は変更しない。

## 1. 結論と正確な範囲

被約実退化因子を持つ rank 四構造に rank 零 plane があり、その全重複度の
abelian-conormal 入口を使う。traceless co-Higgs symbol が非零 global vector field
`V` と一定の非零 normal `sl2` 行列の積であり、trace field `T` と `V` は
線形独立、`[T,V]=0` だが pointwise proportional とする。

この全 source は common-image 四型へ尽くされる。非零 `K` の二型の全射影像は
いずれも二十次元以下。零 `K` の二型では、**全十八 higher 係数を戻したすべての
構造の rank が二以下**であり、rank 四 G source は空である。
したがって本四全源から未知優越成分は出ない。

これは全 G の尽くしでも E の成分帰属でもない。dependent `T=tV`、traceless
symbol 零の scalar source、`Pi0=0` は含めない。normal-nilpotent `b!=0` は
既済の全重複度定理に属し、本定理で再計算しない。

## 2. 全 source の有限化と actual inverse

trace-free lift により `T,V` は `sl3` 内で独立な可換二 plane `A` を張る。
regular 行列を含めば centralizer の Jordan 型により Cartan、mixed、regular-nilpotent
centralizer の三型である。regular 行列を含めず非零半単純元を持てば、その二重固有値
block 内で別の元との一般一次結合が regular になるので矛盾する。従って残りは全元
rank 一 nilpotent である。rank 一行列の二次元線形空間が全体として rank 一なら、
共通 image 又は共通 kernel を持つ。これで common-image、common-kernel の二型を加えた
五型に尽くされる。

標準 ordered basis に対する `det(Euler,T,V)` は順に

| 型 | determinant |
|---|---|
| Cartan | `3abc` |
| mixed | `-3b^2 c` |
| regular nilpotent | `-c^3` |
| common kernel | `c^3` |
| common image | `0` |

である。ordered basis の変更はこの determinant を非零定数倍する。従って今回の
pointwise-proportional 全源は common image だけである。この determinant 表も独立に
全係数再構成した。

common image line `L` を直すと `A=Hom(C^3/L,L)`。その ordered independent pair を
`T=(t1 b+t2 c) partial_a`、`V=(v1 b+v2 c) partial_a` と書く。
`b'=t1 b+t2 c, c'=v1 b+v2 c, a'=a` は determinant
`t1 v2-t2 v1 !=0` の actual base change であり、逆は通常の二行二列逆行列である。
これで `T=b partial_a, V=c partial_a` となる。零 determinant は仮定した pair
independence の外であり、勝手に捨てた source chart ではない。

normal 行列の半単純・nilpotent 型を直すと、`x=(a,b,c,z,w)` において

```
Pi0_ss = partial_a wedge ((b+c)z partial_z+(b-c)w partial_w)
         + k zw partial_z wedge partial_w,
Pi0_nil = partial_a wedge (bz partial_z+(bw+cz) partial_w)
          + k z^2 partial_z wedge partial_w.
```

`k!=0` では次の actual linear maps と tensor scalar により `k=1` へ移る。

| 型 | forward map | tensor scalar | inverse |
|---|---|---|---|
| ss | `(a,b/k,c/k,z,w)` | `1/k` | `(a,kb,kc,z,w)` |
| nil | `(a,b,c/k,kz,w)` | `1` | `(a,b,kc,z/k,w)` |

この Laurent 恒等式は全十括弧と両逆写像で直接確認した。唯一の反転は `k` であり、
`k=0` は別の二型として全て保持する。従って標本点への置換ではない。

共通 line の choice は二次元、ordered basis は四次元、normal semisimple line の
choice は二次元、nilpotent line は一次元である。非零 `K` の係数を戻すと、固定 plane
での affine first-symbol 次元は ss 非零/零で `9/8`、nil 非零/零で `8/7` 以下となる。
この dimension count は全 normal/basis chart を戻した上界である。

## 3. 全 Euler 商と全演算子

quadratic tensor の normal weight 一の raw 源は三十六、linear Euler gauge は六。
既定の三十源とこの全六 gauge の合同行列 rank 三十六を独立に確認した。weight 二は
全九方向であり gauge はない。従ってこれは jet の一部でも、任意に選んだ cut でもない。
quadratic tensor の全 higher source は weight 一三十と weight 二九で尽きる。

各 `Pi0` について全 projective Jacobi 三百五十行への linear operators `D1,D2`、
全六 actual shear の `30 x 6` derivative を再構成し、immutable report の全成分と照合した。
全 kernel に全行を作用させた零と、次の非零 maximal minor で両側の rank を証明した。

| whole 型 | `rank D1 / nullity` | D1 minor | `rank D2 / nullity` | D2 minor | shear rank / minor |
|---|---:|---:|---:|---:|---:|
| ss, `K=zw` | `24/6` | `-3` | `9/0` | `-96` | `6/1` |
| ss, `K=0` | `18/12` | `-1296` | `3/6` | `-27` | `4/1` |
| nil, `K=z^2` | `23/7` | `-2` | `9/0` | `-64` | `6/1` |
| nil, `K=0` | `18/12` | `-81` | `3/6` | `-27` | `4/(-2)` |

非零 `K` では `Pi1` は全 `ker D1` に入り、weight 二の式は
`D2 Pi2 = -Jac(Pi1)` である。`D2` が全九列で単射なので、各 `Pi1` に解 `Pi2` は
高々一つしかない。従って nonlinear compatibility や weight 三・四を課す前から
全 higher fiber は六又は七次元以下。plane choice 六と tensor projectivization 一を
戻すと全射影像は、それぞれ `9+6+6-1=20`、`8+7+6-1=20` 以下となる。
ここでは ss のさらなる空性や既知成分認識を使わず、この最弱の上界で十分である。

## 4. `K=0` の全十八係数と full Jacobi

両型での全十二次元 `ker D1` は、厳密に `Pi1=partial_a wedge W1`、

```
W1 = z((p0 b+p1 c) partial_b+(p2 b+p3 c) partial_c)
    +w((p4 b+p5 c) partial_b+(p6 b+p7 c) partial_c)
    +(B zw+C w^2) partial_z+(D z^2+E zw) partial_w.
```

全十二の表示単位方向が canonical 三十源内で独立であり、全 `D1` kernel であることを
直接検証した。前節の rank 十八と合わせて、これは chosen subfamily でなく全 kernel
である。`W1` は `a` に依存しないので `Pi1` の strict Schouten self-bracket は零。
従って weight 二に非斉次項はなく、**全 `Pi2` は全 `ker D2`** に入る。

その全六次元 kernel は

```
Pi2 = partial_a wedge ((r0 z^2+r1 zw+r2 w^2) partial_b
                      +(r3 z^2+r4 zw+r5 w^2) partial_c).
```

全六方向の独立性と全 `D2` kernel 条件も検証済み。従って全 tensor は
`Pi=partial_a wedge (W0+W1+W2)` であり、全 `Wi` は `a` に依存しない。
全十八係数を同時に独立変数として保持して、strict Jacobi 全十式、Euler-contracted
projective 全三百五十係数、全五 Pfaffian が恒等零であることを別に再構成した。
特に higher compatibility を追加して小さい subfamily だけを残した議論ではない。
全源は rank 二以下であり、actual inverse PGL/scalar を戻してもこの性質は変わらない。

旧 report の linear fiber 十八、像次元三十一/三十は rank 二の大きな全源を数えていた。
その `all_passed=false` は linear-only 判定の不成功として不変に保持し、今回の全構造証明で
rank 四 G source の空性を補う。native GB、jet 延長、追加の標本値を用いない。

## 5. 独立証明書と再現

- 独立 code: `code/classification_plane_global_commonimage_review.py`。
- 有限 job: `research/classification/jobs/plane_global_commonimage_review.json`。
- 独立 report: `reports/CLASSIFICATION_PLANE_GLOBAL_COMMONIMAGE_REVIEW_20261107.json`。
- producer report: `reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_COMMONIMAGE_20261107.json`、
  SHA256 `df0f9f0d2b87780d7caa71a3a8cfc66eacee9fa4d5067b64accf02ee77d381dd`。
- immutable embedded producer source SHA256:
  `ff542fd020159e8409b687fc0423f39f2648691b7fb40c1c77372699e3b1fd68`。

mutable current runner/job には bind せず、report 内の source bytes の hash と今回使う凍結
証拠の hash を検証する。`classification.py gate` は実行前 PASS。独立算術は
3.677 秒、peak RSS 17,340 KiB、全四源 PASS。再現コマンド:

```
python scripts/classification.py gate research/classification/jobs/plane_global_commonimage_review.json
python code/classification_plane_global_commonimage_review.py
```

有限算術の確認と全 source 幾何の読解を合わせた結論であり、gate 単独を数学の証明とはしない。
全分類・全 G・E 帰属・新成分数の結論は含まない。

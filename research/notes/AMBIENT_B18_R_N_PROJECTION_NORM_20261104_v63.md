# Actual R0上でのN projection normと、Rに沿うT0の一般長

2026-11-04。R曲線が物理的にNの外であることだけから、R0 projectionがNの有限像の外にあるとは推論しない。
その不足を、既存の二つの有限unit系の次数2係数だけで直接調べる。

## 1. 同じR0への二つの有限map

前packet `AMBIENT_B18_R_ZERO_SLICE_IMAGE_20261104_v63` は、実際のR曲線と共通19切断から
R0九関数のmapを与える。group19未知数のunit19系を、同じraw parameter Tについて次数2まで解く。
全19式の全三係数が零であることを検査し、九readbackを保存する。

独立受理済み `AMBIENT_B18_N_SAME_R0_BINDING_20261104_v63` のCohen8は、その九関数のうち
実際のrows108,109,112,113,117,118,119,27である。この八つへ上のreadbackを入れ、既存56-unit係数系を
同じく次数2まで解く。全56式を検査する。新しい物理source jetやgroup高次展開は構成しない。

actual N moduleの記号は

```
tau^2=nu, rho61=a+z1*tau, v=rho61-a.
```

したがってR0上のN projectionの正確なannihilatorは fN=v^2-nu*z1^2 である。
literal N8/36 moduleのmonic quadratic relationと、一般点での二次coverから従う。
N tailに付くnilcurveのsupportもこのprojection内にある。

## 2. 四次の最初の非零係数を二次入力から得る

コードはnu(0)=v(0)=z1(0)=0およびnu_1=v_1=0を実際に照合する。
従って

```
[T^4] fN(R0(T)) = v_2^2-nu_2*z1_1^2
```

であり、未計算の三次以上のscalar係数は寄与しない。結果と全19/56補正段階を証明書へ保存する。
この係数が非零なので、実際のR0曲線のgeneric pointはSupp_R0(D4)の外である。
単なる物理的非帰属や、異なるbaseでの古いnorm数値を流用した結論ではない。

## 3. 一般長2の移送

既受理のQ1=Q2=Q3=Jinfとdegree4安定化より、全positive degree piece D_nはD4と同じR0-moduleである。
従って上のR0 generic pointへ局所化すれば全positive pieceが零となり、全relative graded DそのものがD0となる。
さらにD4へのmapも零となるので、同じ局所化でT0=D0である。

実際の全DはR\\N上でordinary doubleであり、Rに沿う一般長は2である。
全positive pieceを消す上記局所化により、この長さはD0とT0へも移る。
nil方向のtorus characterを推測して固定環へ長さを移したのではない。

T0のflatness、R全近傍での二重曲線表示、T0の全curve support、finite socle、元Dの深さは未決。
新しい点・prime・GB・source jet・次数3以上のscalar計算は行わない。

Code `code/ambient_B18_R_N_projection_norm_20261104_v63.py`。
Report/folder `AMBIENT_B18_R_N_PROJECTION_NORM_20261104_v63`。

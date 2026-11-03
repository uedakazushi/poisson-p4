# actual CA/CB odd二次欠損とν非包含：独立検証

`AMBIENT_B18_CA_CB_ODD_QUADRATIC_DEFECT_20261103_v63` を受理する。generic B18で個々の被約actual CA/CB像を \(A=R\oplus It\subset B=R[t]/(t^2-\nu)\) とすると、両族とも **\(\nu\notin I\)**、従ってtracefree \(t^3\notin A\) である。Nで証明した \(I=I+(\nu)\) はCA/CBには成立しない。

## 独立な全source計算

raw CA/CBの物理式を別実装でcanonical化し、20個のfiber parameterを持つ全sourceを構成した。CBでは最初に非線形項 \(-qf\) を含む全sourceを作り、その後で全20変数に関して微分した。この順序により、variable fiber列とflag列の全mixed項、二次項を保持する。

三変数critical tangent planeの係数も、新しい証明書から引用せず、既受理の全critical idealの七つのmonic消去式の一次部分とu,v,wの定義から10×10の線形系を解いて再構成した。全150出力、中心P0、fiber20列、全flag9列からなる150×29 Jacobianが保存値と一致する。

選択した28 actual base行に対して、\(T=T_0+T_1+T_2\) のinverseを直接

\[
R_0=T_0^{-1},\quad R_1=-R_0T_1R_0,\quad
R_2=R_0T_1R_0T_1R_0-R_0T_2R_0
\]

で作り、次数2までの積が単位行列となることを確認した。CBの \(T_2,b_2\) を零とする近似はしていない。全150のSchur微分の次数0、1は零で、二次項は保存した全150成分に一致した。選択28 base行とbranch row61も二次まで零である。producerの一回の新規再実行による全証明書のbyte一致も確認した。

## actual odd二次項への移行

上のvertical vector field Dは28 actual base座標を保ち、\(Dt(0)=1\) である。branch微分と \(\nu=t^2\) から、critical tangent plane上でtは次数3以上となる。全odd一次行列AをこのplaneのRへのrank3微分に掛けると零であり、各b_iのsource一次制限も零である。従って全Schur二次項はactual odd係数の二次制限を与える。even部分やtrace shiftの寄与を無視した推測ではない。

core座標pとflag座標pのtracefree基底の違いは、Rの単位によるtの変更で、その中心値は1である。full R firstjetが零のresidualの二次class、およびこのcritical tangentへのchosen liftの二次制限は、その変更に影響されない。

## 全150生成元に対する欠けるclass

actual R29の全150 odd一次行列はrank7である。七つの基底行への全一次従属係数を再計算し、全b_iからその定数結合を引いた。得られた**full R firstjetが零**の全150 residualについて、二次制限の係数を全て計算した。

両族とも、その空間は正確に \(\operatorname{span}(WU,WV)\) である。UV、U²、V²、W²の全係数が零で、row46、95のWU/WV minorが非零である。UVを加えた係数空間のrankが3となることも確認した。

一方、一次項が正確に \(d\nu\) となるactual I内の定数結合 \(i_\nu\) の二次制限は、CAで \(-UV/75\)、CBで \(-16UV/9\) である。これら全係数は保存証拠と一致した。

仮に \(\nu\in I\) なら、\(i_\nu-\nu\in I\) はfull R firstjetが零である。そのplane二次制限は、\(\nu=t^2\) が次数6以上であるため、上記の非零UV項を持つ。

しかし任意のfirstjet零のIの元 \(\sum f_i(R)b_i\) について、plane二次項に寄与するのは \(f_i(0)\) だけである。b_iのsource一次項が零なので非定数係数との積は次数3以上になる。定数係数はAの一次関係を満たし、従って全150 residualの定数spanで尽くされる。その二次制限はWU/WVのspanに入り、UVにはならない。矛盾であり、両族で \(\nu\notin I\) が全次数のidealの非包含として従う。

この議論は、cubic以上の未知生成元を暗黙に省略していない。それらも同じ有限生成idealの係数付き和に含まれ、上の二次制限を変えられないため、高次計算でUV欠損を修復することはできない。

## 保持する未解決範囲

全Iの生成元、cubic class、minimal生成数の上界、depthやsubintegrality、actual Dの完成型は決定していない。\((s_6,\nu+cUV,WU,WV,U^3)\) のような候補を採用しない。

R19\Nでは各CA/CBは個々に滑らかで、normalizer旗は一つである。CAの旗とCBの旗が異なることを、同じ成分のnormalizationの二preimageと扱うconductor議論は使用しない。今回の結論はactual odd係数の独立な有限計算から従う。

独立codeは `code/ambient_B18_CA_CB_odd_defect_review_20261103_v63.py`、receiptは `reports/AMBIENT_B18_CA_CB_ODD_QUADRATIC_DEFECT_20261103_v63_REVIEW.json` である。

```sh
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python \
  code/ambient_B18_CA_CB_odd_defect_review_20261103_v63.py \
  --replay work/B18-CA-CB-odd2-producer-replay \
  --output work/B18-CA-CB-odd2-independent-review.json
```

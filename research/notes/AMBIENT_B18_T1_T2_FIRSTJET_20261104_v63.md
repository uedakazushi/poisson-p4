# Q1/Q2のactual honest9入力と、共にrank43の一次関係

2026-11-04。新しいsource jetやHadamard係数を作らず、凍結済みQ3の正確な入力から
Q1/Q2をliteral column subsetとして取り出した。

H=R0+m9*tauの36Koszul関係を保持する。
full degree分解では

```
Q1 = lower degree0*H + degree1*C0   （1513列）,
Q2 = lower degree0,1*H + degree2*C0 （1926列）。
```

両方の各normalized pairはQ3の同じpairに一致し、そのhonest10 liftも同じcached membership
規則である。全pairのpointer対応を確認し、全liftを再利用した。

唯一のunit関係graph52（weight1、H multiplier1）は、Q1の0-based index1243、
Q2の1353に実際に入っている。全36Koszulを先頭に置けば同じ10/3 pivotで
free generator5を消せる。したがってQ3と同じ9×10射影を正確に適用できる。

```
Q1 honest9 presentation relation columns: 1548,
Q2 honest9 presentation relation columns: 1961。
```

この射影と全残余列を保存した。元のrelationのうち単位同値で消した一列以外を落としていない。
全一次係数はFULL9 FIRSTJETの同じcached scalar-function入力に対応する列を取り出したもの。
四つのg商の未知微分の感度は両方の全列で零であり、一次matrixはexactである。
E0の全一次部分も零。

正確な有理数Gaussにより、両方の一次関係spanは **rank43** と分かった。
独立column indicesと全90係数を保存した。
これは43本の一次部分についての結果であり、残る二次classを含めた最小関係空間全体ではない。
従って、このpacket単独でT1/T2の零・非零、cyclic性、全Dのdepthは主張しない。
Jinfの最小44本の定理と二次functionalを合わせる操作は後続の別gateである。

Source `code/ambient_B18_T1_T2_firstjet_20261104_v63.py`、
`reports/AMBIENT_B18_T1_T2_FIRSTJET_20261104_v63/Q1-Q2-input-and-firstjets.json`。
新しいprefix、IFT、source jet、GB、点は無し。一回2.087878秒。
全2069Q3一次入力を再計算せず、そのactual subset bindingを証明した。

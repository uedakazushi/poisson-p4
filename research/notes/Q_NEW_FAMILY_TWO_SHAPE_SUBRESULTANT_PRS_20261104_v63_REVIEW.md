# 二 shape 結果式の独立 Brown PRS identity audit

2026-11-04。実施するのは一回の120秒/2GB/32MiB制限付き検証である。
原著の FLINT resultant の再実行、別 prime、点、gcd、source26 は使わない。
係数環は別の reversed LEX 順 `QQ[Z2,Z1,Z0,R,r]` とし、s 多項式を短い係数配列で扱う。
二入力の284/1076全係数を受理済み三制約から組み直し、整数化200/400を確認する。

`A=E2` のs次数5、`B=E1` のs次数3に対し、予想される通常次数列なら

```text
C2=-prem(A,B),
C1=prem(B,C2)/LC(B)^3,
C0=prem(C2,C1)/LC(C2)^2.
```

これは Brown subresultant PRS である。全除算を係数多項式上の完全除算として扱い、
積へ戻す恒等式も確認する。各 pseudo-division は
`LC(B)^d*A=Q_DAG*B+R` を与える。Qは各段の
`Q_next=LC(B)*Q+LC(R)*s^(degR-degB)` というDAGで保持し、巨大なBezout係数を展開しない。
初めに小さい monic/nonmonic 二例について、別の8×8 Sylvester整数行列式と比較する。
原著入力順は `(E1,E2)` なので、次数3×5の奇数符号により最終比較は
`raw_Res(E1,E2)=-C0` となる。

途中で予想次数2→1→0を外れる、時間・メモリ・証明書容量の上限に達する、または完全除算が
失敗した場合は、その段階までを保存し停止する。別方式へ自動移行しない。
**最終全係数比較が完了するまで、47,114係数の結果式 identity を独立受理したとはしない。**
完了状態と最後の段階はreport/native-frontierに記録する。

degree1 subresultant `a*s+b` が保存され、aが単元である開集合では `s=-b/a` と回収できる。
この追加開条件を消したり、rの一意回収・actual NEW target体のindex1・実定義方程式を
主張したりしない。原著のnamed boundary stripは独立の多項式分解であり、本検証の
結果式 identity が未完了なら、その分解だけから identity を補わない。

```sh
source /workspace/cloud-setup/activate.sh
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/q_new_family_two_shape_prs_review_20261104_v63.py \
  --output work/v62-review/two-shape-independent-prs --work work/v62-review/two-shape-independent-prs-build
```

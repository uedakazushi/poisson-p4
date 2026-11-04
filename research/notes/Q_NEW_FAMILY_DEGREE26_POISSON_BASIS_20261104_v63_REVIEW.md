# Poisson生成Inv26基底の独立ambient微分・係数検証

2026-11-04。同じ旧三ベクトルの形式平面を使うが、微分は**制限前の全16個の
ambient beta座標**で取る。C,Dの全二次係数からunit16方向を先に微分し、
その一次式へ固定平面を代入する。二つのplane変数だけの微分を代用しない。

独立Python疎整数辞書engineで五trace4・三trace6の全8×16gradientを作る。
各wordのcyclic product ruleとEuler恒等式を照合し、元のOmega逆を全16方向
に収縮する。全10のdegree14 bracket、3本F4Inv14と15bracketからなる全18の
degree18候補を再構築する。theta/Hの微分や追加candidateは使わない。

degree26の57列は、Inv14三基底×受理済み全Inv12九基底の27列と、primitive
Inv18六基底×trace Inv8五基底の30列。Inv12の値のみは自分の既凍結の独立
全基底係数証明を再利用し、今回の著者係数行列をコピーしない。57のsemantic
定義と全係数を独立に検査する。

minorは3、9、45角で7584、26750、7998 mod32003。さらに元の16重みから
完全対称冪DPとC2 Weylの全8項を別に計算し、Inv14/18/26の次元3/9/45を
確認する。p単元なQQ係数の非零minorとこの上界より、選んだ45本は全Inv26の
QQ基底である。character upperをformal planeのrankから推定していない。

Sp4はOmegaを保つため、ambient定数Poisson括弧は不変式を不変式へ送る。
degreeはdegf+degg−2であり、全候補の普遍性と次数が保証される。N12とraw
trace(CD²CD²)の差はF4Inv8に属するため、その置換はF4Inv14成分を加える
だけ。F4はmoment Casimirなので不変式とのPoisson括弧は零である。

この結果は普遍基底だけ。actual NEW点、source26の展開・関係式検査、degree26
方程式の存在は含まない。既存のNEW次数下界26は変わらない。

```sh
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/q_new_family_degree26_basis_review_20261104_v63.py --output work/degree26-basis-independent
```

60秒・2GB。著者producer/nativeは再実行せず、一回の独立係数計算で確認する。

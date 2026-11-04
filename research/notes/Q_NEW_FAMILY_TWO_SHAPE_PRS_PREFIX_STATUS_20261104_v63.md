# 独立 PRS の正の prefix と最終段階の停止

2026-11-04。`Q_NEW_FAMILY_TWO_SHAPE_SUBRESULTANT_PRS_20261104_v63_REVIEW`
の一回のジョブは14.243秒で停止した。**all_passed=false** を保持する。
最終 pseudo-remainder 計算中、FLINT が1,075,879,976 bytesの追加確保に失敗してabortした。
プロセスには2,000,000,000 bytesのaddress-space制限を課していた。
代替計算や元resultantの再実行は行っていない。

完了した有限証明は次の通りである。

- 全284/1076入力係数を受理済み局所化式から独立に構成し、整数化200/400を確認。
- monic/nonmonic二つの小さいfixtureで独立8×8 Sylvester整数行列式と比較。
  行列式は549と208997で、原著の入力順の結果式が `−C0` となる符号も一致。
- `C2=−prem(E2,E1)` はs次数2、全係数項数23,788で保存済み。
- `C1=prem(E1,C2)/LC(E1)^3` はs次数1、全係数項数54,980。
  除算は係数ごとのexact divisionと再乗算で確認し、全係数を保存した。
- 二つのpseudo-divisionとcofactorのDAG recipe、leading factor、入力・出力hashを保持。

従ってC1は、表示したleading factorが単元の開集合上で、実際の両入力idealに入る
一次式である。別の独立packet `Q_NEW_FAMILY_TWO_SHAPE_COPRIME_DEGREES_20261104_v63_REVIEW`
は分母の互いに素性とmap degree3/5を証明し、generic incidence上でsを有理回収できることを
示す。その開集合では保存されたC1から `s=−C1_constant/C1_linear` を使える。

一方、最終C0は未保存・未検証である。元の47,114係数raw resultantとの全係数一致、
r射影19次の実成分の認定、actual target体の指数1、NEW定義式を主張しない。
以後の既知因子除去やmemory-reduced normは、新しい明示的な上限付き許可の対象とする。

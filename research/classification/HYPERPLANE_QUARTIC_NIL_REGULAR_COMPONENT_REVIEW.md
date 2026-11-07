# regular nilpotent 四次源：全周囲の一法線支持の独立確認

2026-11-07。課題 `E_reduced_nonnormal`。ここでは一つの22次元全空間成分の存在と、
証明点での被約成分の唯一性を確認する。既知89成分との相違、成分数の増加、
他のcofactor型の被覆は別の義務であり、本レビューから昇格しない。

## 1. 実際の点と独立算術

全nil4・translation零の線形源に
`n1=1,n8=-1/4,n12=-1/3`、他のnを零として代入し、hで斉次化して36倍した点を使う。
座標は `(h,x,y,z,w)`、hはCasimirで、非零六括弧は

    Pxy=36hw+6xz−8y², Pxz=6wx−10yz, Pxw=−10wy−z²,
    Pyz=6wy−9z², Pyw=−12wz, Pzw=−12w².

実Pfaffianは `9((z²−2yw)²−48hw³)` で非零である。
四次因子はhについての原始一次式であり、w³と `(z²−2yw)²` は互いに素なので既約である。
したがってこの点は単純超平面と既約四次因子を持ち、一般rank四の実Poisson点である。

独立実装は `code/classification_quartic_regular_component_review.py`。
以前の独立nil22 reviewerの標準ライブラリFraction・疎多項式演算のみを再利用し、
producer、legacy producer helper、SymPy、NumPyをimportしない。
全源の六式に上記パラメータを別途代入して点を再構成した。
全十strict Jacobiが零であることを確認し、その後に全150係数を自由な周囲変数として扱った。

以下を3.34秒、最大RSS 23,312 KiBで確認した。

* 全350個の射影Jacobi係数と、全150列の微分Dを独立に構成し、`rank D=102`。
* 実Euler鉛直25列、実GL5作用25列、scale列を別構成し、その全像Mのrankは47。
  鉛直だけのrankは24であり、射影軌道次元は `47−24−1=22`。
* 公開された一法線Rについて `DR=0`、`rank(M,R)=48`。
  `dim ker D=150−102=48`なので、これは全kernelを張る。
* 全二次ベクトル `Q=Jac(R)` を再構成し、`rank(D,Q)=102`。
* 独立な102次minorで `DT=−Q` をFraction消去し、公開Tと全150係数が一致。
  選択行だけでなく、全350行で `DT+Q=0`。
* `B(R,T)=Jac(R+T)−Jac(R)−Jac(T)` の全350係数を再構成。
* 全kernelを張るMとRの全列に対し `B(R,ker D)⊂im D`。
  全mixed列をDへ加えたrankは102のまま、さらに `B(R,T)` を加えると103となる。

最後の非零性はrank宣言だけに依存しない。
102次minorの102行に元方程式の144番行を加え、独立な左汎関数λを構成した。
λは全D列と全kernel mixed列を消し、厳密に

    λ B(R,T) = 9108855/2 ≠ 0

となる。103次minorと非零の定数determinantも証明書に保存した。
このλは公開producerの余核汎関数を読み込まず、再構成Dから独立に解いたものである。

## 2. 真正な線形切片と三次障害

Dの確認済み102列に対応する標準座標方向の空間をUとする。
102次minorが非零なのでDはU上単射である。
全kernelは `M⊕CR` だから

    C150 = M ⊕ (U ⊕ CR), dim M=47, dim(U⊕CR)=103.

従ってS=`U⊕CR`は、実際の47本の線形cutで定義される切片である。
単にkernelを群方向で割ったベクトル空間と、真の非線形Poisson支持を混同しない。

選択した102方程式をUについて形式陰関数消去し、残る座標をtとする。
全解の初めは `P+tR+t²T+O(t³)` である。
TはUに属する唯一の二次補正で、全350式で二次係数が消えることを確認した。
三次係数への新しい補正の寄与はDの像だからλで消え、
λが取り出す残式は `(9108855/2)t³+O(t⁴)` となる。
係数体は標数零なのでこれは非零の三次式に単元を掛けたもので、真の切片の被約支持は原点だけである。

さらに、仮に二次補正を別のTへ変えても、その差は全kernelに属する。
その三次係数の差 `B(R,ker D)` は全てDの像と独立確認したので、λの値は変わらない。
この完全な補正曖昧さの確認は、固定切片での議論に加えた強い検証である。
正確な非被約環やnil指数の同定は本レビューの結論に加えない。

## 3. 全周囲成分と唯一性

47本の真正cut後の支持が零次元なので、Krull次元不等式から元の150係数支持の
証明点での局所次元は47以下である。
実軌道・鉛直・scaleの像が47次元なので下界も47で、一致する。
Euler鉛直24方向と非零scale一方向を引くと、全射影Poisson空間の局所次元は22である。
実PGL5軌道も22次元だから、その閉包は全空間の22次元既約成分となる。

唯一性を次元一致だけから推測しない。
連結なPGL5は、全Poisson空間の有限個の既約成分をそれぞれ保つ。
従って証明点を通るどの全空間成分も、同点の22次元軌道閉包を含む。
同点での局所次元上界22と、既約成分の極大性から、それらは全て同じ軌道閉包に一致する。
よって証明点での被約ambient成分は唯一である。

この唯一性を使えば、同点を含む既約な全sourceの像閉包をこの成分へ入れることができる。
ただし全sourceの正しいパラメータ化・全境界の保持は独立の幾何的義務である。
その入口は `HYPERPLANE_QUARTIC_NIL_REGULAR_SUPPORT.md`。
本レビューは個々のsource全体が一軌道であるという、より強い主張を使わない。

## 4. 再現と結論の範囲

```sh
python3 scripts/classification.py gate research/classification/jobs/hyperplane_quartic_regular_component_review.json
.venv/bin/python code/classification_quartic_regular_component_review.py --output reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NIL_REGULAR_COMPONENT_REVIEW_20261107.json
```

出力先には既存証明書を上書きしない別名を指定する。
gateは契約とhashを確認しただけで、上記の算術実行・数学的含意とは別である。
全データ、実minor、独立λとpairingは上記REVIEW reportに保存した。
既知89成分との相違・新成分数・残る全分類はrootの別監査まで未決とする。

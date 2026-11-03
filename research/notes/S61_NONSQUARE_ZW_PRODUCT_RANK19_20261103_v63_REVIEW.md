# 元の ZW product19・generic quartic146 の独立監査

2026-11-03。対象は `S61_NONSQUARE_ZW_PRODUCT_RANK19_20261103_v63`。以下の全域下界19と generic 等号146を区別して受理する。

元の六二次式の保存された全係数から q0,q1,q3,q5 を再構成し、Z=3v0−v1+v6=0、W=v3=0、heavy=0 のみを代入する。五つの独立変数は v0,v1,v2,v4,v7 であり、ell をさらに零にはしない。全四二次式の係数を著者証明書と比較し、全 aC-open で

`q1=L²/8100`, `q3=L(eL−400a*x0)/8100`, `L_x1=−40C`

を確認する。非零 universal linear syzygy があるため rank20 は成立しない。

q0 は x1 を含まず、L=0 後の四変数対称行列式を、係数多項式の24項 determinant から独立計算すると `625a²/1024` になる。従って代数閉体上でも q0 modL は既約である。L=x0=0 の q0,q5 pair の全12係数を元式から計算し、5a の単位 pivot 後の5成分と保存された30項の scalar row を比較する。各 column の正規化比は a,C と非零QQ定数だけの unit である。保存された cofactor を自前 Fraction 辞書で掛け合わせ、literal sum=1 を確認する。117保存項には零多項式を表す一つの零項が含まれ、非零項数は116である。この数え方の相違は恒等式に影響せず、著者の frozen bytes は変更しない。著者の native checker、pivot producer、GB は実行しない。

よって pair は全 aC-open で独立であり、q0,q5 modL は互いに素。四式を unit 行変更で `(Q0,Q5,L²,L*x0)` とすると、任意の linear syzygy を modL、次に mod(L,x0) へ落とすことで Q0,Q5 の係数は零になり、最後は L,x0 の一つの Koszul syzygy だけが残る。L,x0 は L_x1 の可逆性により独立である。従って20本の一次倍は全 aC-open で正確に rank19を持つ。d、J の逆元を用いない。

generic CI は別の exact certificate で確認する。元804項から a=e=d=1 の qfour modell を再構成し、保存四変数切断で degree5 の56×80 Macaulay 行列を作る。指定56角 minor の非零QQ determinant を独自 Fraction Gauss 消去で再確認する。これは一般点で barJ が六変数CI4であることを示し、全点CIとは言わない。

元六式に対する `I=J+ell*(Z,W)` の unit 行変更と同一座標の認識は、独立受理済み RAW24→product41 bridge の source/係数 hash を依存として固定する。一般点で barJ がCI4なら J もCI4で ell は S/J の nonzerodivisor。従って

`(I:ell)=J+(Z,W)`、

`0→S/(J,Z,W)(−1)→S/I→S/(J,ell)→0`

が成立する。六変数CI4の次数4は48、五変数ZW商の次数3は35−19=16なので、一般点で HF4(S/I)=64、七変数 light 四次の像は210−64=146である。

元の係数行列は QQ[a,e,d] 上の多項式行列である。一般点で rank146であることから全147角 minor が恒等零となり、全ての多項式特殊化で上界146が従う。物理 chart の意味は従来の aFBC-open に保つ。全点の下界146を標本や抽象自由性から推論しない。

heavy×六二次式は light 四次と別の列にある。元六二次式の6角QQ minor を同じ非零点で独立計算して generic rank6を確認したので、重み四 quadratic-products は generic rank152、全点上界152である。これは147/153を146/152へ訂正する結果であり、四次全射・長65・ambient除去を示さない。

再現（新しい出力先が必要）:

```sh
source /workspace/cloud-setup/activate.sh
python code/s61_nonsquare_ZW_rank19_review_20261103_v63.py \
  --output work/marked-square-v63/ZW19-independent-new
```

全 frozen archive/input hash を照合し、60秒・2GB以内の検証だけを行う。新しい体、GB、有限体からの昇格、cofactor の探索は行わない。

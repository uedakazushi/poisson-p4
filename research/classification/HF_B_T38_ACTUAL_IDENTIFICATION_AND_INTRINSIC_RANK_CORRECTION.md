# HF-B全源はT38に含まれる：実PGL同定とrankの訂正

2026-11-08、PRIMARY E。以下は同じ固定anchorと完全有限graphに対する実恒等式であり、新しい点やprofileの試験ではない。独立レビューは別義務である。

## 実座標変換

HF-Bは `HF_Km2_E0_WALL_2*p + 1` の全Poisson A10源である。既証明の像の上界23と固定anchorの実軌道次元23から、PGL5×A10の像閉包はそのanchorの軌道閉包に等しい。anchorを通る完全graphはt次数三であり、全350式の全係数が恒等零である。真の102次元sliceと形式陰関数定理の一意性から、このgraphは唯一の滑らかな射影24次元成分を定める。

保存した整数化t=3eのe=1を用い、元座標をy=(a,b,c,u,v)、新座標をX=(x0,x1,x2,x3,x4)=Myとする。

```
M = [0  1   0  -1/6   0]
    [6  0   9   -1    3]
    [6  0  -9   -1   -3]
    [0  0  -1    1    1]
    [0  0   1    1   -1]
```

行列式は−288であり、保存した有理逆行列は左右の逆を満たす。完成済みT38源のパラメータを

```
(p0,t0,t1,t2,t3)=(24,1,−1,12,−12)
```

とすると、tensor scalarは1で、実raw恒等式は

```
M_* Π_HF-B(e=1) = Π_T38 + E_X ∧ M_*R,
R=(12c+4v, 0, 4a, −2c−2v, 3a−4u).
```

全150二次係数、全10のE∧残差、および全Euler-Pfaffian五次式をQQ上で検証した。符号を含む式は

```
D_T38(X)=det(M) D_HF-B(e=1)(M^{-1}X).
```

これは射影Poisson構造の実PGL同定である。T38の閉PGL不変成分はこの点を含むので、唯一の24次元局所成分はT38である。graphの非零eは正規化トーラスにより同じ実PGL軌道に属し、e=0のanchorも閉T38に入る。最後に前述の全A10像閉包の等式から、HF-Bの全源がT38に入る。全A10の各点について新しい軌道正規形を仮定していない。

## 以前のrank推論を撤回する

`HF_ANCHORS_SINGULAR_SUPPORT_AND_NO_SECOND_PLANE.md` のBに関する「S′の一般点でraw Π_uv=−2u²が非零なのでPoisson rank二」という一文は誤りである。射影rankをrawの一成分から判定できない。この推論を明示的に撤回する。同noteの古いbytesは保存し、今回の訂正を外に置く。

実際、元anchor上のtwisted cubic cone

```
(a,c,u,v)=(−2s*t²,t³,3s³,3s²t), b自由
```

へE∧Πの全10三次式を引き戻すと全て恒等零となる。従ってS′は射影Poissonの零surfaceであり、一般rankも零である。Dのsingular supportがS∪S′であり第二planeがないという記述はそのまま成立する。rawrankの誤推論は全D101、実群rank48、真のslice、トーラス重み、完全graph、または全A10閉像の証明には使っていない。

## 検証範囲

`CLASSIFICATION_PLANE_M2_HF_B_T38_ACTUAL_IDENTITY_20261108.json` に、Mと逆、元と先の全bracket、明示Euler場、全150零残差、全10零残差、両Pfaffianと変換式、cone上の全10零式、および実行source全文を保存した。0.548秒、QQ完全算術。jobのgateは実行前PASS。共有runner/helperを編集していない。

HF-Aの親の同定は未解決である。今回のT38同定からHF-Aに関する結論を推測しない。

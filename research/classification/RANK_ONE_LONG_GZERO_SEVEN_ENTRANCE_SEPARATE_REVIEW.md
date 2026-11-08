# g零七whole源の25次元入口：独立pure review

2026-11-08、PRIMARY E_multiple_hyperplane。**独立pure PASS。** producer `RANK_ONE_LONG_GZERO_SEVEN_25_IMAGE_ENTRANCE.md`（SHA256 721f956e545e96bdcb48944cfd7d7596b5d2c8240b2d9958a4ff03471171f091）の全源・全次元・全prior／相互分離を、元六括弧とaffine作用から別に導出した。

新point、M1、Jacobi operator、ambient matrix、数値rank、kernel、jet又はCASを実行していない。以下は**25次元閉族の入口**の受理であり、七族の全空間成分性は未証明。

## 1. 全g零源と有限表

元(19)でg零なら全条件は(U+1)H=0、Fは任意二次式六係数。PxyからF、Pxz/PxwからU、Pyz/PywからdHを戻せる。定数差も(U+1)c=c=0により消えるのでmapはinjectiveである。Hを選んだ点又は一つのFへ固定していない。

producerのB1–B6を九つのmonomial重みiλ+jμ=−1へ代入した。各表で二つの表示monomialが全固有空間を尽くし、他のdegree三以下のmonomialは固有値−1を持たない。λμ(λ−μ)は全六行で非零。qの中心を実affine translationで消せ、その逆も全sourceへ戻る。A/B零は元五monomial型へ戻す。

V1はU=−z∂z、H=z(Aw²+Bw+C)の全eigenspaceを持つ。A又はdiscriminant零を一般点として流用しない。leading零はlinear／constant pへ、double-rootはz*(linear)²へ戻り、元M1/M11/M21に被覆される。

## 2. 七つの正確な12→22→25次元

B1–B6はU1の相異固有値conjugacy orbit二、translation二、全固有空間二、全F六で12。eigenline二つのcoefficient unit開も既約であり、有限frameの選択を余分な連続parameterとして数えない。

V1はU1のrank一・trace−1 conjugacy orbit二、image内translation一でU三。各Uの全H-eigenspace三、全F六で12。Uのzero-eigenlineには任意の中心があるため、Uだけから中心二を取る数え方は不可である。producer最終版はこの点を修正し、直接3+3+6又はHの二rootの中点をmarkした4+2+6だけを使う。どちらも余分なfiberを持たない。

g零のStabX作用を別に導出した。r′=r+Δ(q)、q′=Bq+tなら

```
U′=B_*U, H′=det(B)H,
F′=F−U(Δ)−(1+d)Δ
```

（全q-pullbackを戻す）である。全七sourceはこの作用及び逆を保つ。normalized dilationはU′=U,H′=λH,F′=λ²Fであり、sourceはcone安定。

δP=Xのequivarianceと係数逆からAff4作用写像の全fiberはStabX十次元のcosetで尽くされる。従って全affine像は20+12−10=22。全表示開でDの唯一の重複hyperplaneはhなので、H-choice四とscalar一を戻し、Hを忘れる一般fiberは一点。七つとも射影像閉包が正確に25次元となる。

## 3. actual factorと相互識別を元式から確認

Pf P=−Hを斉次化した全七Dはproducerの表と一致する。

* B1はh³とrank三quadric、他六族は最大factor重複度二。Q rank三はA/B unitから直接従う。
* B2/B5のcubicはそれぞれz／hについてprimitiveな一次式で、gcd条件から絶対既約。全singular supportはそれぞれplane{h=w=0}／plane{z=w=0}。前者だけがmarked Hに含まれる。単なるSing(C)∩H非空という常に起こり得る条件は使わない。
* B3/B4/B6はunique double-H、unique simple-L、unique irreducible rank三Qを持ち、直接制限すると(rank Q|H,rank Q|L)=(2,1),(1,2),(1,1)。重複度がHとLを区別するので交換は出来ない。
* V1は三つのsimple因子z,w−ρ1h,w−ρ2hを持ち、ρ1≠ρ2。simple三formのspanは三。SCのsimple三formは中心を通るbinary pencilのspan二に入るため相違する。

全marked factor／singular-plane／restriction-rankはprojective multiplicationとGrassmannianのproper incidenceへ戻せる。同次元閉既約族が等しければ、各dense constructible像の非空開が交わることも用い、genericに相反する上記不変量から七者の不一致を得る。

## 4. 全prior及び既新像からの分離

旧prior129の全projection、rank二六、rank四dimension二十四以下百十四という入口を保持する。これらはrank又はdim25下界で除外。CA/CB二十八の全Dはall-linearなのでB1–B6を含まない。V1にはCAの三乗factorもCBの二つの二乗factorも無く、proper掛け算像の外である。

旧nil31三25は同じall-linear／反復factor条件で排除。positiveNEW25のgeneric Qはrank五であり、B1の全rank≤3 quadric閉条件と一致せず、他六族のgeneric Dにはtriple-Hが無い。T01/T39/T40のgeneric被約開と七像の全multiple-H閉条件も両立しない。同次元包含が等号を強制する場面と、より大きいCA/CBの全閉条件による排除を混同していない。

非定数g新25のgeneric Q rank五に対し、B3/B4/B6は全Q rank≤3、V1はall-linear、B1はtriple-H、B2/B5はirreducible residual cubic。従っていずれも同一25像ではない。

全七Dはessential≤3のclosed determinantal条件を満たす。定数g二共鳴のgenericD essential五も独立に確認した。Cのpartial matrixにおけるproducerの二determinantは、表示された順序で対角積−4g0³t²及び−12g0³tcとなる。D=h²Cの零方向はmod hでまずh方向を失い、次にCの零方向となるのでDもessential五。従って新二共鳴との同次元一致は無い。

SCはall-linearでB1–B6と異なり、V1とは§3のsimple-form spanで異なる。M1/M11の全像はall-linearかつ四乗／三乗factorを持つので、同じgeneric比較で七像全体をそれらのsourceへ帰属させることも出来ない。

## 5. 受理した範囲

全七sourceの正確な25次元、各全zero/unit routing、全prior129／非定数g／新三源／M1/M11との分離及び七者相互分離を受理する。旧sourceの全134／350方程式を再実行したというreceiptではない。

fullambient D100、新成分数、generic軌道及び全E完了は本reviewの結論ではない。次の有限componenthood判定が必要な状態を保持する。producer本文・manifest及び既数値証明書のbytesは変更していない。

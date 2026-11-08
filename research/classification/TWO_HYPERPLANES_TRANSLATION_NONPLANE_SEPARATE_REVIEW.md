# 二超平面・非平面rank二translation二全源の独立MPS被覆

2026-11-08、PRIMARY E、`/root/log_manuscript`。
全source・flag/orbit閉包の幾何と、essential四／Perazzoの二つの実SS変形を
別作者として独立確認した。**二全源のMPS40被覆をPASS**とする。
新成分、ambient局所一意性、七法線の支持又は全1+1+3被覆は主張しない。

算術は既独立Fraction ringを使う `classification_two_quadrics_snc_review.py`
のisolated `--translation-nonplane` 一回だけで再構成した。
producer、S61 engine又はCASをimportせず、producer埋込sourceはhashの検査にだけ使用した。
実測0.232639146秒、最大RSS 23932 KiB。新point、anchor、normal jet、GBは無い。

## 1. 二つの必要十分な全source

固定cofactorはV=x∂u+h∂v、h-Casimir代表を用いる。
`TWO_HYPERPLANES_CUBIC_TRANSLATION_RANK_TWO_NONPLANE_SOURCES.md`
にある全21係数の五括弧を、表示式から独立に組み立てた。
元strict Jacobiの非自明項は(u,v,w)一つだけであり、その全base三次係数を
旧nilpotent reportの15式とmultisetで完全一致させた。
従って保存された残式の一部だけを必要条件として用いたのではない。

A=n11、B=n18とすると三次因子の非平面部分は
`(Ah-Bx)(hu-xv)` である。A非零ではactual w'=w/AでA=1、
B/A非零ではさらにx'=(B/A)x、u'=(B/A)uでA=B=1へ戻る。
いずれも元{x,u}=x²、{x,v}=hxを保ち、反転する尺度以外のparameterを捨てない。

A=0,B非零ではactual swap(h,x,u,v,w)→(x,h,v,u,w)を行う。
swap後のh行はh(0,0,x,h,0)であり、
`E wedge(0,0,x,h,0)`を引きtensorを−1倍するとh-Casimir代表と同じVへ戻る。
三次の非平面部分は `(Bh-Ax)(hu-xv)` となるのでA,Bが交換される。
これは実PGL、Euler代表変更、非零scalarだけを用いた全chartの変換である。
A=B=0の平面源は別補集合として残す。

両branchで全15式をliteralに代入し、producerの各三角stepの
元式・pivot・係数・解をFractionで再導出した。
L=hでは11step、L=h−xでは12stepである。
L=hだけn17を反転するが、n17零ならF|h0が零でF全体がhで割れるため、
既約三次を仮定する当該全源ではn17非零が必要である。
ほかのstepは非零数値係数だけで割っている。
最終表示を全Jacobiへ戻して逆方向も恒等的に確認した。

従ってfixed exact VのL=h源は非零Aの一尺度と
`A7×Gm(n17)`の八自由係数を持ち、次元9。
L=h−x源は非零A,Bの二尺度と全A7の七自由係数を持ち、次元9。
全normalizerと逆、det1、全parameterのactualD=hxF、変換後のcofactorを確認した。
この次元9はchosen pointのtangent値ではなく全supportの有限chartから得る。

## 2. 射影flagの全13次元とfaithful fibre

射影Πとordered二factor H1,H2から
`X=Π sharp dlog(H2/H1)` の射影linear field[N]を取る。
rank四symplectic openではdlog(H2/H1)が非零なのでXも非零。
trace零代表Nはrank2、N²=0であり、
二hyperplane covectorsはP(im N*)=P1の相異なる二点である。

I=imNの選択はGr(2,5)の6次元、I⊂K=kerNの選択は
K/I∈P(V/I)の2次元、V/K→Iのprojective isomorphismは3次元。
ゆえに[N]全体は11次元、ordered二点はさらに2次元で合計13。
これはGL5が推移的に作用するflag空間であり、無断のgeneric flags限定は無い。
順序swapは有限なので次元を変えない。

一つのflagをh,x,[V]へ送ると、projective Πのscaleを一意に決めてVをexactに固定できる。
h-Casimir代表にEuler ambiguityも残らない。
実際E wedge Lのh行が零ならhLj=xj Lh。
linear Lとhの整除からLh=ch、Lj=cxjとなりE wedge L自体が零である。
従って正規化後のprojective fibreは§1の全21係数sourceそのもので、高々9次元。

各全非平面sourceのPGL飽和像閉包Sに対し、これで **dim S<=13+9=22**。
固定Nのcentralizerを選んだsliceのdimだけから全像へ飛ぶ推論ではない。

## 3. 同じ既宣言点の実22orbitとwhole閉包の等号

`TWO_HYPERPLANES_TRANSLATION_NONPLANE_CANONICAL_MODELS.md` の
同じ二つのkappa=1原source点を全source表示から評価した。
全25Euler列E wedge(xj∂i)、全25actual GL列
`AΠ+Π A^t-xj∂iΠ`、scalar列を全150係数で独立に再構成した。
保存済みfullFの先頭51列と全て一致し、Fraction消去で

```
rank Euler=24,
rank(Euler,actual GL,scalar)=47
```

を確認した。さらに八又は七sourceの全微分も独立に組み立て、
fullF全行列と一致しrank47を得た。
よってprojective actual orbitは47−24−1=22次元である。
ambient D96又は七法線の再計算はこの最弱の帰属証明には不要なので実行していない。

normalized八源A7×Gm及び七源A7は既約、PGL5も連結既約である。
従って各Sは既約で、上の22orbit Oを含み、§2によりdim S<=22。
Oの閉包も既約22次元なので厳密に **S=closure(O)** となる。
この等号が一点の陽性を全sourceへ運ぶ橋である。
一点の陰性からSの全境界の陰性を導くものではない。

## 4. 真の二SS親とSNC成分定理

座標(H,X,U,V,W)のH行を零にする。
Perazzo用は全A,L,a∈A3に対し

```
XU=2XU, XV=XV, XW=-XW,
UV=aUV, UW=-LH V-(a+2)UW, VW=AHX+aVW;
D=HX[2AHXU+LH V²+2(a+1)UVW].
```

essential四用は全A,b,a∈A3に対し

```
XU=3XU, XV=XV, XW=-XW,
UV=aUV, UW=-bV²-(a+3)UW, VW=AHX+aVW;
D=HX[3AHXU+bV³+3(a+1)UVW].
```

左辺XYは括弧{X,Y}を表す。
両全親のstrict10、元projective350、全actual Pfaffianを独立に展開して零／表示Dと一致させた。
PerazzoではA=L=0、essential四ではA=b=0、a=1で
それぞれ4HXUVW、6HXUVWの真のSNC rank四点を持つ。
全A3の像閉包は既約なので、その全像を含む一つのambient成分はこのSNC点にも会う。
既受理の全SNC-meeting component theoremにより、そのambient成分はMPS40の一つ。
従って各親の全像がMPS40有限閉和に入る。
SNC一点が既知成分に属することだけから、任意の族を逆向きに所属させた論証ではない。

## 5. 二つのactual arc、全kappaと全pole相殺

Perazzo arcはinternal座標(h,x,u,v,w)で

```
U=u+h/ep, V=v+x/ep, W=w+(h-x)/ep,
a=-2+2kap*ep, A=ep^-2-2kap/ep, L=-2/ep+4kap,
tensor scale=ep,
G(h,x,u,v,w)=(h,-2x,-2v,u/2,u/2-v+w).
```

essential四arcは

```
U=u+h/ep, V=v+x/ep, W=w+h/ep,
a=-3+6kap*ep, A=2ep^-2-6kap/ep, b=-6ep,
tensor scale=ep,
G(h,x,u,v,w)=(h,x,v,-w,u/6+w/2).
```

全shearのdet1、両Gのdet−2及び−1/6、両向き全逆をFractionで確認した。
元の三つの親parameterをこれらのliteral式に代入し、parameter積を独立化していない。
得た全150係数でep負冪が完全に消え、全ep,kapについてstrict10と元350が零となる。
全actual quinticもproducer表示と係数ごと一致した。

ep零端点は、両方とも既宣言の同じPi_kapの全150係数へ厳密に一致する。
essential四のfixed prior変換v_new=v_old−2wも逆へ戻して、元kappa=1 source点との全150一致を確認した。
kapは一度も反転せず、kap零もこの正則曲線に直接含む。
端点Dはそれぞれ

```
Perazzo: hx[(h-x)(hu-xv)+hxw],
essential四: hx(h²u-hxv+x³)
```

で非零、従って全kapでrank四である。

ep非零では実invertible PGLと非零scalarを施した真の親の点だからMPS40有限閉和内、
ep零では全tensorが正則に延びるので端点も同じ閉和内。
特に同じkap=1 canonical orbit Oがこの閉和に入り、PGL安定性と閉性から
§3の全S=closure(O)も入る。
これでessential四の完全八parameter源、Perazzoの完全七parameter源、
A=0/B非零swap chart、全非零normalizer及び各sourceの全境界を被覆する。

## 6. 保存と残余

新独立reportは
`CLASSIFICATION_TWO_HYPERPLANES_TRANSLATION_NONPLANE_SEPARATE_REVIEW_20261108.json`。
jobは `two_hyperplanes_translation_nonplane_separate_review.json`、
入口manifestは `two_hyperplanes_translation_nonplane_separate_review_manifest_20261108.json`。
旧checker全FunctionDef ASTと旧source／二review契約bytesは
`two_quadrics_snc_pre_translation_review_source_contracts_20261108.json`へexactに保存した。
旧producer、旧report、旧契約及びmainには変更を加えていない。

平面A=B=0、他nilpotent/cofactor、残three plane1+1+3源、Ga9linear-w、
stable源、全G/E及び分類主定理はこの新結果から完了とはしない。
二22orbit閉包がambient成分であるとも、親の局所一意性があるとも主張しない。

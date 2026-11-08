# Jordan(2,1)、(L,M)=(h,x)全六源：actual SS confluenceの独立確認

2026-11-08、PRIMARY E、`/root/log_manuscript`。**全六源のMPS40被覆PASS**。
producer `TWO_HYPERPLANES_J21_H_X_SS_CONFLUENCE.md` と全packetを別作者として確認した。
算術は既独立Fraction checkerのisolated `--j21-hx` 一回で行い、
producer、S61 engine又はCASをimportしていない。実測0.072861秒、最大RSS 22328 KiB。
新point、anchor、normal、jet又はGBは使っていない。新成分数は零である。

## 1. 全六源をliteral inverseで二parameter模型へ戻す

`TWO_HYPERPLANES_CUBIC_SUPPORT.md` §6及びaccepted source dataの
`j21_h_x_irreducible` は、元の全45係数から得た必要十分な六source入口として使う。
今回その前提となる全14自由係数のbracketsを表示式から組み直し、
strict Jacobiの全三次係数を保存された12式とmultisetで完全一致させた。
accepted六source bindingを代入し、全strict10／元350を恒等零と確認した。
この確認を、元14係数の他枝を新たに尽くした証明とは呼ばない。

全六parameterはn0,n5,n7,n8,n10,n12、既約三次stratumはn5非零である。
A=n12、B=n12*n8+n7、C=n0+n10*n12としてactual変換

```
U=u-Aw-Bx-Ch, V=v-Ax, W=w-n10h-n8x
```

を独立に組み立てた。det1、全逆、両方向の全150係数、
cofactor V0=v∂u+h∂v+x∂wのexact保存を確認した。
source parameterは一つも反転していない。
s=n5、k=n8とすると結果は全域で

```
xu=xv, xv=hx, xw=x²,
uv=v²/2-hu, uw=s h²+kxv-vw, vw=khx-hw,
D=hx[-s h³+x(v²/2-hu)].
```

左辺ijは括弧{i,j}を表す。
s非零ではcubicのu係数-hxと残部xv²/2-s h³は互いに素なので、primitive linear polynomialとして絶対既約。
s零のreducible境界も下の同じcurveには残している。
四つの余分な係数を指定点へ置いたsliceではなく、全inverseを持つ全源正規化である。
optional amplitude scaling又はdense orbitの推論はこの被覆に使わない。

## 2. 全三parameterの真の親とSNC橋

H行零、自由(A,b,a)∈A3に対して

```
XU=3XU, XV=XV, XW=-XW,
UV=aUV, UW=-bV²-aUW, VW=A H²+aVW,
D=HX[3A H²U+bV³+3aUVW]
```

の親を独立に再構成した。全三parameterを記号のまま残し、
strict10、元projective350、全actual Euler-Pfaffianを確認した。
A=b=0,a=1でD=3HXUVWは真のSNC rank四点である。

親A3の全像閉包は既約であり、その全体を含む一つのambient成分はこのSNC点にも会う。
既受理の全SNC-meeting component theoremにより、その成分はMPS40の一つ。
したがって全三parameter親がMPS40の有限閉和へ入る。
SNC一点の既知帰属だけから任意の曲線へ所属を逆向きに推論していない。

## 3. 全s,kを保持する一つのactual曲線

e非零で

```
H=h, X=x,
U=w+x/(3e), V=u+v/(2e)+h/(2e²), W=v-h/e,
a=1-3ke, A=a/(2e³), b=-8s e⁴,
tensor scale=e
```

を真の親へ代入する。
座標mapのdet1と全逆を確認し、特にA=a/(2e³)の積関係を維持した。
e以外を反転せず、得られた全150係数のe,s,k冪は全て非負となった。
全e,s,kのstrict10、元350、全actual Pfaffian及び親五次のactual pullbackを確認した。
e零の六括弧は§1模型へ係数ごと厳密に一致する。
従ってnormalizerの逆により元全六parameter源の端点一致も全150方向で戻る。
k零、s零及び全sourceの他係数零をこの実curve自体に含む。

e非零は真の親の点のinvertible PGLと非零tensor scaleなのでMPS40閉和内にあり、
e零ではcurve全体が正則に延びるため同じ閉和内。
全six source、全inverse正規化、PGL飽和及びその閉包が被覆される。
ambient接次元、局所一意性、normal obstruction又は新成分性は不要である。

## 4. 保存と補集合

reportは `CLASSIFICATION_TWO_HYPERPLANES_J21_H_X_SS_SEPARATE_REVIEW_20261108.json`、
jobは `two_hyperplanes_j21_h_x_ss_separate_review.json`、
manifestは `inputs/two_hyperplanes_j21_h_x_ss_separate_review_manifest_20261108.json`。
旧checker全FunctionDef ASTと三current契約のexact bytesは
`inputs/two_quadrics_snc_pre_j21_hx_review_source_contracts_20261108.json`へ保存した。
旧source、producer、report、契約及びmainを改変していない。

j21(h+x,h) p零全五源、j21(0,h+x)のu依存六源及びn6零補集合、
他nilpotent型、残る三plane1+1+3源、全E/G及び分類主定理は別に残す。

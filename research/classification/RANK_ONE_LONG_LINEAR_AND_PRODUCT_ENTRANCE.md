# rank一long・g零M1/M11：二つの全25次元像と全priorからの分離

2026-11-08、PRIMARY E_multiple_hyperplane。**純全source入口。新点もambient計算も未実行、全空間成分性は未証明。** CA/CBのmarked factor境界を正しく保持し、係数型の違いだけで28次元親を排除しない。

## 1. 二つのwhole sourceと正確な次元

`RANK_ONE_TRANSLATION_WHOLE_FINITE_EXHAUSTION.md` のg零全源を使い、常にFは全二次式六係数である。

* **M1**：Hは非定数affine linear式。Hの三係数に対してU(H)=−Hを課すと、Uの六係数にはrank三の線形条件となる。H=zにactual q-affine正規化した全解はU=−z∂z+v(z,w)∂w、全affine v三係数。従ってH3+U3+F6=12のinjective既約sourceである。一般Dは4H∞+L、H∞が唯一の重複hyperplane。
* **M11**：Hは相異なる非平行二affine linear factorの積。係数空間は二本のaffine line各二、非零乗数一で既約五次元。H=tzwにactual正規化した全Uはλz∂z+(−1−λ)w∂w、全λ一係数。従ってH5+U1+F6=12。一般Dは3H∞+L1+L2で三つは相異なり独立、H∞が唯一の重複hyperplane。

係数からH,U,Fを一意復元する全逆、Stab(X)十次元と全action-fiber十次元、scalar stabilityは既源入口と同じ。両sourceはStab(X)を保ち、一般点ではmarked H∞をDから一意に復元できる。従って両PGL飽和射影像閉包Z1,Z11は**正確に25次元**となる。pointの軌道rank又はβ像を次元下界として使わない。

## 2. marked Poisson hyperplaneのmodular行列はregularに依存する

任意の二次lift Pからcanonical divergence-free lift

```
Πbar=P−E∧div(P)/5
```

を取る。これはEuler gaugeの選択によらず、Πbarの係数はprojective tensorに線形に依存する。Poisson hyperplane ℓ=0に対してY=Πbar#dℓをcolumn-sharpで取る。

projective tangencyはℓ零上でYがEulerに比例することを意味する。次数と線形座標のgcd一から

```
Y=ℓA−A(ℓ)E
```

と書け、Aはlinear vector field。δΠbar=0からdiv Y=0なので

```
0=ℓ*tr(A)−5A(ℓ).
```

従ってY=ℓB、B=A−tr(A)E/5はtrace零、B(ℓ)=0。ℓ-Casimir liftはΠbar+E∧Bであり、そのdivergenceは

```
N_ℓ=5B=5Πbar#dℓ/ℓ.                         (1)
```

N_ℓはtrace-free homogeneous linear fieldとして一意である。ℓのscaleでも変わらない。ℓの一係数をunitにしたdual-projective chartでは、多項式のℓによるexact divisionはregularなので、N_ℓはmarked tensor／hyperplaneにregularに依存する。極を仮定した無界なframe normalizationを使わない。

従って「ℓ^m divides actual D、ℓ Poisson、rank N_ℓ≤1」はmarked factor incidence上の閉条件である。ℓはproperなdual P4で動くため、その存在条件もrank四Poisson開集合で閉じている。D非零の乗法写像によるfactor incidenceを用い、Pf零へ落ちた点から不正に因子を選ばない。

## 3. CAとCBの全marked係数

### CAのgeneric triple-H

既CAのh-Casimir普遍源ではPxy=f(x,h)、Pzw=g0h²、残るyq項は任意三次Hamiltonianの微分である。全modular場はf_x∂yであり、そのhomogeneous行列はrank一以下。一般Dは3H+H1+H2なので、triple-Hのmarked incidenceとrank N_H≤1がCA全rank四閉包へ延長する。

M1/M11ではℓ³|Dとなる唯一のhyperplaneは元H∞。そのmarked modular場は

```
N_h=h∂x+x∂y,
```

でrank二、Jordan型(3,1,1)。従って両全像はCAに含まれない。

### CBの二つのgeneric double-H

`MULTIPLE_HYPERPLANE_LOW_RANK_MODULAR.md` §4.2 の全rawCB源を使い、座標(h,X,y,z,w)で

```
f=X(h−2X)/(6a), a≠0,
Pxy=f, Pzw=μf,
Pyz=−f_X*z/2−H3_w,
Pyw=−f_X*w/2+H3_z, Pxz=Pxw=0.
```

任意三次H3、全磁気方向を保持したままdiv P=2f_X∂yである。二つのdouble factorに対する**row** cofactorは

```
L1=X: Y1=(h−2X)∂y/(6a),
L2=h−2X: Y2=−2X∂y/(6a).
```

P−E∧Yiがそのfactor-Casimir liftなので、全divergenceは

```
N_L1=(−3h+2X)∂y/(6a),
N_L2=(2h+2X)∂y/(6a).
```

どちらのhomogeneous行列もrank一以下で、H3にもμにも依存しない。(1)のcolumn-sharp規約とrowの符号を混同しない。

両double-labelを残したproper factor incidence上でこの閉条件を延長できる。M1の一般D=4h*Lでは両generic double-labelがhへcoalesceし、simple-labelはLへ行く。唯一ℓ²|Dのhyperplaneもhなので、境界でrank一条件が失われることはない。実際は一つのdouble-labelだけを保っても同じ帰結となる。N_h rank二との矛盾でCBもZ1を含まない。

M11の一般D=3h*L1*L2は、そもそもCBのproper掛け算像D=H*L²*M²の外。marked-rankを追加しなくても排除できる。

## 4. 他の既知成分と二像の相違

rank二全成分は閉rank≤2、rank四dimension24以下はdim25下界で除外する。残る旧・新25次元成分への包含は、同次元なら一致を意味する。

Z1の全rank四部分はproper閉条件D=H^4*Lに入る。既知25成分にはgeneric Dがこの4,1型のものは無いので一致しない。Z11は全D all-linearかつtriple-Hの閉条件を持ち、MPS25のgeneric被約D、positiveNEW25のgeneric smooth-Q、非定数g成分やSCのgeneric最大重複度二、定数g二共鳴のgeneric既約cubicと異なる。nil31 DDHのgeneric2,2,1型にはtriple-Hが無く、half/thirdの全Dは3,2型のproper閉条件なので、Z11のgeneric3,1,1型はその外。

Z1とZ11も、前者の全Dが4,1閉条件内で、後者の一般Dがその外であるため相異なる。これは両25像の分離であり、より大きな未知成分への包含を排除したものではない。

## 5. 未解決の正確な範囲

両sourceの全coefficient、全λ／affine v、全Fと全zero/unit境界は上の全源に保持される。しかし全ambient成分性又は未知親はまだ決めていない。新しい固定点、rank profile、M1又はjetは未実行。現在認証したのは、正確な25次元像二つと全既知閉成分からの分離までである。

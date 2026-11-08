# Ga9 linear-w：全T38と新A5成分を親候補から除く

2026-11-08。PRIMARY E。新しいpure geometry。CAS、point、matrix、minor、profile又はjet無し。
Ga9 c零linear-wの全十一源の飽和閉包は、既受理の一つのactual23軌道Oの閉包である。
以下はその**密軌道が全T38及び新A5 smooth-hyperflex23成分に入らない**ことを示す。
全E帰属は未解決であり、新成分性又は全境界点の非所属は主張しない。

## 1. 正しいcanonical tensorとその全projective zero support

`NIL211_SHORT_TRANSLATION_REVIEW.md` と保存magnetic normalformにより、
全源閉包の密軌道を次の整数scale四の代表で調べれば十分である。
座標(h,x,y,z,w)、h-row零、非零entryは

    Pxy=4h²−y², Pxz=−hx+4w²+2yz, Pxw=16hy−4wy,
    Pyz=hy, Pzw=−4hw+4y²。

PGLによるprojective rankをraw一成分で判定しない。
intrinsic zero schemeは全十のcubic E∧Pで定義される。
そのうちh²y及びh(4h²−y²)から、任意のzero pointでh零が従う。
h零上では、残る四つのcubicは

    −y(4w²+3yz), 3wy²,
    4xy²+6wyz+4w³, 4y³。

したがってy零、続いてw零。
逆にh=y=w零では全十cubicが消える。ゆえに厳密に

    Zero(P)_red = {h=y=w=0} ≅ P1。                 (1)

nilpotentな厚みが存在してもKrull dimensionは支持と同じ一である。
全23軌道上でこのdimensionは一に保たれる。

## 2. completed T38の全一般源にはrank零surfaceがある

入力は受理済み `MPS_SMOOTHING_COMPLETION_20261103_v62/T38.json` の
**全correctionを含む**五parameter族である。
p=p0として、p t2 t3非零の稠密open上で

    q1=(p/2)x1x3+t2x4²,
    q2=−(p/2)x2x4+t3x3²,
    q0=p x1x2+(4t2t3/p)x3x4

と置く。completed sourceの全非零entryは正確に

    P01=−2t0q2/p, P02=2t1q1/p,
    P12=q0, P13=q1, P24=q2。

q_iを同時に消す三次曲線は、四座標(x1,x2,x3,x4)において

    [s:t] ↦ [−2t2 s³/p : 2t3 t³/p : s t² : s²t]

というrational normal cubicである。
x0自由のconeはdegree三のintegral surfaceであり、その全域でrawP、従ってE∧Pは零。
t0,t1を制限しておらず、p,t2,t3の零補集合へは次の閉性で戻す。

## 3. zero surfaceの存在は全成分閉包で保たれる

全projective bivector空間上で、universal zero scheme

    Z = {([P],q) : (E∧P)(q)=0} ⊂ P(H0(Λ²TP4))×P4

を考える。これは全十cubicの閉零集合で、第二因子P4がproperなのでZからparameter空間への射影はproper。
従ってfiber dimension上半連続性により

    { [P] : dim Zero(P)≥2 }

は閉集合である。
同じ結論は、cubic-cone surfaceの固定Hilbert polynomialに対するproper Hilbert incidenceからも得られる。
zero schemeのflatness又はboundary conductorの型一定は仮定しない。

MPS全族定理により、§2の稠密parameter openとPGL飽和の閉包は全T38である。
上のclosed条件はPGL不変なので、**T38全体がdim Zero≥2を満たす**。
(1)よりGa9 linear-wの密23軌道は全T38に入らない。
これは一つのnormal-cone chart又は一つのsparse degenerationの陰性ではない。

この閉条件は他のtargetにも適用できるが、それぞれで全intrinsic zero supportのdimensionを
別に証明する必要がある。rawentryだけのrank判定へ戻さない。

## 4. 新A5 hyperflex成分は同dimensionの異なるdense orbitである

Ga9のcanonical actual divisorはD=hF、

    F=4h³w−8h²y²+y⁴。

Fはintrinsic plane quotient上で唯一E6 unibranch cuspを持ち、H={h=0}はそのcusp tangentである。
この全source／絶対既約性／23軌道閉包は先行独立定理を使う。

一方、`E6_A5_TWO_NORMAL_COMPONENTS.md` とその独立全二法線証明が受理した
A5 smooth-hyperflex component C_A5は、全六magneticを含むactual23軌道OAの閉包である。
そのplane quarticは

    F_A=h²y²−2hw²y−16hy³+w⁴

で、A5の二つのsmooth branchesを持つ。Hは別のsmooth pointでhyperflexとなる。
Ga9のE6一枝とはanalytic branch numberから既に異なる。

仮にGa9の密点がC_A5に入れば、PGL不変性によりclosure(O)⊂C_A5。
両者は既約dimension23なので等号となる。
代数群の軌道は局所閉であり、closureにdenseな軌道はその中のopen集合である。
二つの23軌道O,OAが同じ既約閉包内のopen集合なら交わり、従って同じ軌道でなければならない。
しかしDは唯一のlinear factorHとirreducible quarticFを持つため、実PGL同値はこのfactor pairを保つ。
essential rank三のgradient kernelからintrinsic P2 quotientも定まり、そこでquartic cuspの枝数も保たれる。
E6とA5の違いに矛盾する。

したがって密Ga9軌道は**新A5の全ambient component**にも属さない。
quarticのsingularity adjacencyだけから、同dimension componentの全sourceを親と呼ぶことは出来ない。
逆にこの非所属を、Ga9の全magnetic-zero境界がC_A5と交わらないという主張へ拡張もしない。

## 5. 残す義務

本稿で新たに除いたのは全T38及び全C_A5への密軌道包含である。
既済の全T39/T40 graph陰性、103 marked-pair排除等は各凍結証明のscopeを維持する。
T36、T37又は他のhigher-dimensional親のactual recognitionは未解決のままであり、
二次支持P∪Lをactual曲線又は新成分と呼ばない。
旧二つのnonclosing局所route、固定N twist、monomial pencil、jet列を再開していない。

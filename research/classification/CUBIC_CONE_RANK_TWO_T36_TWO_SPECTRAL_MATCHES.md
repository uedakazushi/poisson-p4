# T36比較における二つのspectral scalarの明示

2026-11-08。凍結済み `CUBIC_CONE_RANK_TWO_ONE_ORBIT_FRONTIER.md` §4の
正負二対応を、符号と変数順を省かず記録する純証明追補。
新operator／point計算は無い。

T36のaffine quadratic modular fieldは
`(3p0/5)*diag(−3,3,1,−1)`。
我々の整数代表の同じfieldは
`5*diag(−3,−1,1,3)`、変数順は(x0,x1,x2,x3)。
T36 tensorに掛かる任意非零scalarをalphaとすると、固有値集合の一致は

`alpha*p0=+25/3` 又は `alpha*p0=−25/3`

だけである。実数の大小で順序を選んだのでなく、集合
{−3,−1,1,3}を非零scalarで保つものが正確に±1だからである。
相異なる四固有値なので、その対応を満たすlinear mapは次の変数順の各線のscaleのみ。

* 正の場合は(y1,y2,y3,y4)=(x0,x3,x2,x1)。
  我々のjetの二entryは
  `Q12=+5*y1*y2+10*y3*y4, Q14=+5*y1*y4`。
* 負の場合は(y1,y2,y3,y4)=(x3,x0,x1,x2)。
  その二entryは
  `Q12=−5*y1*y2−10*y3*y4, Q14=−5*y1*y4`。

どちらもy1*y2とy1*y4のlogarithmic係数比は1。
T36のliteral completed sourceでは全t0,…,t3について、同じ比が
`(alpha*p0)/(alpha*p0/5)=5`。
固有線をyi→di*yiでscaleしても、bracketのdi*djとmonomialのdi*djが相殺し、
logarithmic係数は変わらない。全tensor scalarも比では相殺する。
従って二つの可能なspectral scalarは両方とも除外される。

前noteの一意order≥2零点の証明により、もしactual PGL同値が存在すれば
T36のdisplayed vertexは必ず我々のc-vertexへ送られる。
vertex固定PGLのquadratic jetへの作用は上のGL4のみで、vertex shearは高次だけを変える。
よってこの比較は一点同士の恣意的なmarked比較に依存しない。
結論はT36 generic dense orbitとの非同値であり、higher-dimensional成分の
specialization又は未証明componenthoodをここから主張しない。

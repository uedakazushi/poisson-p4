# B18/W全曲線の固定Darboux変換と境界の有理葉

2026-11-03。generic B18代表P0と、三者のcritical W曲線の共通物理式を使う。
同次座標(h,x,y,z,w)で

    P_W = ∂y∧(x²∂h+hx∂z+(h²+z²)∂w) + W*x²∂z∧∂w.

全CA/CB sourceの保存された150係数と一致する。これは既存Nの局所完成型の記号tとは別の、実Poisson二ベクトルのbirational geometryである。

## 全Wを保つ固定polynomial Darboux写像

x=1の開集合で座標を(h,y,z,w)と書く。非零bracketは {h,y}=−1、{y,z}=h、{y,w}=h²+z²、{z,w}=Wである。

    H=h, Y=y, Z=z−h²/2,
    T=w−h*z²+(2/3)h³*z−h³/3−(2/15)h⁵

と置くと、全Wに対し {H,Y}=−1、{Z,T}=W、他は零となる。逆式は

    h=H, y=Y, z=Z+H²/2,
    w=T+H*Z²+(1+Z)*H³/3+H⁵/20.

Wを反転しない固定多項式自己同型で、Jacobianは1。射影斉次化は次数5、逆も次数5である。全五同次座標の両合成はx²⁴倍の恒等写像となり、共通因子がないのでCremona bidegreeは正確に(5,5)である。

## 軌道と葉

g_r=diag(1,1,r,1,1)について、全150係数で r⁻¹ g_{r*}P_W=P_(W/r)。従ってW≠0は単一のprojective PGL軌道である。generic Wで独立受理された局所完成型があれば、この恒等式で全W≠0へ移せる。ここで軌道式自体の証明は完成型の仮定を必要としない。

W≠0では開集合は定数symplectic A4の単一葉である。W=0ではZ,Tが二つの多項式Casimirで、全affine葉はH,Yを座標とするA2となる。従ってgeneric B18点P0の葉もこの開集合で有理である。B18の全特殊軌道にこの結論を拡張しない。

Z=c,T=dの葉の閉包はy方向をvertexとするconeである。その底曲線は、P1の[S:A]から

    [h:x:z:w]=[A*S⁴:S⁵:c*S⁵+A²*S³/2:
       d*S⁵+c²*A*S⁴+(1+c)*A³*S²/3+A⁵/20]

で与えられる。S≠0ではx成分が非零、S=0ではw成分が非零なのでbasepoint-freeである。h/x=A/Sが有理逆parameterとなるため曲線はrational quintic、leaf closureは次数5の有理coneとなる。c,dの全有限値でこの議論が成立する。

コードは全bracket、両逆式、同次両合成、全W軌道式をQQ上で検査し、元の物理sourceへのbindingも保存する。新しいPoisson成分や成分数の増加、全global分類は主張しない。

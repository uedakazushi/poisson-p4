# v=0全uに対するT39/T40 graphの一変数厳密判定

2026-11-04。元のactual P(u,1,0)を、固定されたaffine原点で正常化する。
新しい数値β/source点は取らず、全計算をQQ[u]で行う。一回60秒/2GB/32MiB。

固定値はuに依存せず、Pfaffianは−16/81。constant Darboux frameは

    R=[[1,0,0,0],[0,0,0,9/8],[0,1,0,−1/32],[0,0,−9/2,0]],
    det R=−81/16.

値をJへ移した後に、凍結された歴史的S61.normalizeのselected-raw perspective gaugeを使う。
primitive Lee gaugeと取り違えない。perspective responseの逆はconstantで、全61座標に
u依存の分母はない。全150 reconstruction、元134方程式、projective JacobiをQQ[u]で確認する。

T39とT40には、既受理の**全A16上の閉じた多項式45 graph** Q51とQ52を使う。
discriminant openや選択group-minor openへの所属判定ではない。
actual q(u)−Q51/52(beta(u))の45式ずつを全て保存し、QQ[u]のmonic gcdを計算する。
従ってu≠0の各componentへの所属は、保存されたgcdからuの因子だけを除いた多項式が
零となることと厳密に同値である。定数gcdなら全u≠0の非所属が従う。

具体的gcdと全90多項式はreport／graph-residuals.jsonへ保存する。
u=0は今回のu≠0境界分類の対象とせず、gcdのu-adic位数だけを記録する。
T01の全点排除とT39/T40の一般排除は、先に凍結した別の閉性／orbit次元packetを使う。
この計算はgeneric非所属を点の標本から推測するものではない。

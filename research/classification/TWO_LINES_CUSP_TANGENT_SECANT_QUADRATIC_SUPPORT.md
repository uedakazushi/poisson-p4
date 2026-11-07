# 尖点接線と尖点割線：未計算だった全六法線二次写像

2026-11-07、PRIMARY `E_reduced_nonnormal`。新成分又は実変形は主張しない。

## 既済の計算と今回の有限対象

`TWO_LINES_PLANE_CUBIC_SUPPORT.md` の平面十三配置で残る三配置のうち、
`node_flex_tangent_through_flex` と `cusp_flex_secant_through_flex` は既に
D101/source48の全一法線二次が零、固定三次も零である。
また1+2+2の残二配置は既に全二法線二次・三次が非孤立である。
これらを再実行せず、旧runnerの法線1/2の選択条件から除かれていた一配置だけを扱った。

対象は `C=Y²Z−X³` と二直線 `Y,Y−X` の全六磁気源
`cusp_cusp_tangent_cusp_secant` である。
旧reportで固定されたr=(1,2,3,5,7,11)を変更していない。
全150係数・全350射影Jacobiの微分はD96、Euler/GL/scalar像は47、
全磁気source微分を加えた実像Fは48で、真の余法線は六次元である。

## 全二次写像と非閉鎖の証明書

既存runnerを無変更で `--normals --normal-dimension 6` として一回だけ実行した。
全六法線のSym²に対応する21列を、全350行で計算し、
Dの96独立列・96独立行の非零有理blockで余核へ射影した。
その全行空間は七次元で、Dの全150列を消す七つの有理行汎関数と、
全法線、全七二次式、exact grevlex基底26本をreportに保存した。

全七二次式は t2=t3=t4=t5=0 へ制限すると恒等零である。
従って二次支持は少なくとも二次元平面を含み、原点支持ではない。
これは全二次写像の結論であり、選んだ二方向だけの検査ではない。
一方、各literal法線軸のJacobiは非零であり、この平面を真のPoisson曲線又は
全成分の次元下界へ読み替えることはできない。

## 真切片から全成分へ必要だった矢印

F48、法線N6、Dの96独立標準方向U96は全150の直和となる。
従ってC=span(N,U)はFに横断する真正102次元の線形切片である。
96implicit変数を非零D blockで解けば、残六変数の二次leading idealが上記全二次idealとなる。
もしこれが原点支持なら切片次元0、Krull不等式でraw次元48以下となり、
実全source像48の下界と合わせてsource閉包がprojective23成分だと示せた。
Fは実群47にsource方向を追加したものなので、48次元群との形式的積は仮定しない。
今回は前提の原点支持が成立しないため、この矢印から成分性を結論しない。

## 再現と停止範囲

有限契約は `jobs/two_lines_cusp_tangent_secant_quadratic.json`。
元のsource/全係数を保持したexact reportは
`reports/CLASSIFICATION_TWO_LINES_CUSP_TANGENT_SECANT_QUADRATIC_20261107.json`。
gate PASS後2.393秒で完了した。旧producer・共有Ga9runner・旧証明書は変更していない。

これで残る平面1+1+3三配置、1+2+2二配置は全て、既点での全二次支持が
成分性に十分でないことが明確になった。高次延長・法線の部分probe・別anchorは行わない。
全源のactual変形又は退化因子の幾何へ経路を変更する。

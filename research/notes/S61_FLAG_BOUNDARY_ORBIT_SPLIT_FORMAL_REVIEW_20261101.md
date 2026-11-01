# 疎な旗境界の二方向分離：独立レビュー

2026年11月1日。
[主note](S61_FLAG_BOUNDARY_ORBIT_SPLIT_FORMAL_20261101.md)の固定β命題を承認する。
指定点で `0≠t∈m⁴\m⁵`、`t²=0` である。
相対ambientや十八次元境界全体の除去はこの結果からは結論しない。

## 独立計算と再現の区別

独立コードは表示された有限行列から群法則、symplectic恒等式、指定βの保持を直接確認する。
最高双ベクトルを有限行列でpullbackし、射影chartを再正規化した後で微分する。
二接方向のv0,v1 minorは4であり、他の七接成分は零である。
この最後の座標比較には主producerと同じraw局所frame utilityを使う。

元135式のQQ四次jet、observable tの次数判定、QQ Gröbner分解と有限長moduleの
Hilbert関数は主producerを別に再実行し、8.77秒でPASSした。
これらを独立アルゴリズムで再実装したとはしない。

## 局所環の論理

有限群は実際にβを保ち、全135式の同変性がある。
指定点で群接方向と元の自由座標による線形切断v0=v1=0が横断的なので、
群作用写像はambient chartでétaleとなる。
従って固定β完成環は切断環S'と二形式変数の積である。
tは射影再正規化の単元倍となるので、次数・非零性・平方零性は切断で判定できる。

初期部分ideal J1に対し、v2飽和B1は自由階数5、E=B1/J1は有限長111、最高次数6である。
同じ実在旗族の切断曲線C'は滑らかで、b1微分が `−v2/2` であるため、
一般点は既証明の自由5・t=0の旗開に入る。
有限横断商によりS'はR=QQ[[v2]]有限、v2はreductionであり、上界5と下界5が一致する。

追加初期核Kは階数零なのでEに含まれる。
R-torsion TにはS'のm-adic filtrationを**誘導**する。
torsion元を殺すv2級数の最低次項はその初期元を殺すので

    gr(T) ⊂ E/K

となる。右辺の全次数上界6により `T∩m⁷=0` である。
主noteはK=0という未証明の等式を使っていない。
一般点でt=0なのでt∈T、QQ四次observableの非零剰余から `t∈m⁴\m⁵`。
よってt²は `T∩m⁸=0` に入る。

observableを第136方程式としてidealへ追加してはならない。
主コードは元135式から選んだpivotだけで陰関数を消去し、追加行からt座標を読む。
この区別を確認した。

## 再現

```
.venv/bin/python code/s61_flag_boundary_orbit_split_formal_review_20261101.py --output /tmp/flag-split-review.json
.venv/bin/python code/s61_flag_boundary_orbit_split_formal_20261101.py --output /tmp/flag-split-replay.json
```

主証明書、主note、ソースのSHA256は独立レビュー証明書に固定する。

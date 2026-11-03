# 全 R19 は既存 LPB39 に入り、全 Poisson 環はその各点で非 CM である

2026年11月3日。これは全交叉の成分R19を、既存の一般階数二Poisson成分LPB39へ同定する包含である。
新しいPoisson成分を加えるものではない。

## 定理

全閉集合R19は既存LPB39に含まれる。
従ってRの **全ての閉点** における、全P4 Poisson空間の完成局所環はCohen–Macaulayでない。
これは個々の被約R,CA,CBの局所環や、二成分交叉だけの局所環とは別の主張である。

## 1. 既存LPBの定義との一致

`legacy/v14/code/rank2_lower_bounds/lpb_family.py` はLPBを次の通り定義している。
射影中心[e]∈P4と、商C5/Ceの任意のhomogeneous quadratic vector field Vを取り、

    P=e∧V

とする。Vはcenter方向に依存しない。
radial fields ell*E4の四方向はEuler-verticalとなるので、vector-field parameterは40−4=36次元である。
centerの四方向と合わせてaffine40、projective39となる。
この族が全Poisson空間の実際の成分であることは、既存の全Jacobian tangent upper40と
family lower40のcertificateで証明済みである。
本作業ではこの既成功の大きなJacobian計算をやり直さず、元のsourceとcertificateをhashして入力にした。

Rの八係数族は元から

    P=∂y∧V(h,x,z,w)

であり、Vはyに依存しない二次場である。
座標yを旧LPBのe4へ移す置換(0,1,4,2,3)で、元の **全八係数** の双ベクトルは
LPBの表示に生の行列成分でもcanonical150係数でも厳密に一致する。
LPBは閉でPGL5不変なので、この八係数族の全PGL飽和閉包Rも含まれる。
一般点だけの所属ではなく、特殊化と境界を含む全Rの包含である。

## 2. 全Poisson完成局所環の非CM性

任意の閉点p∈Rを取る。pは既知の全Poisson成分CA28,CB28,LPB39の全てに入る。
これらが全Poisson空間の実際のirreducible componentsであることを既存の成分証明から用いる。
従ってpでの全Poisson局所環には、商の次元が28となるminimal primeと、39となるminimal primeが存在する。

有限型C代数はexcellentであり、完備化後もこれらのcomponent dimensionsは保持される。
もし完成局所環がCohen–Macaulayなら、complete local ringのunmixednessにより全minimal componentsは
同じ次元を持たなければならない。28と39の差はこれに矛盾する。
従って全完成局所環は非CMである。

この論証は全Poisson成分の尽くしを必要としない。
他に未知成分があっても、既知の異なる次元の二つのminimal componentsだけで結論できる。
R自身がPoisson空間のminimal componentだと扱ったものでもない。

## 証明書と依存

* `code/ambient_R19_LPB_containment_20261103_v63.py`
* `reports/AMBIENT_R19_LPB_CONTAINMENT_20261103_v63/report.json`
* 既存定義: `legacy/v14/code/rank2_lower_bounds/lpb_family.py`
* 既存full-component certificate: `legacy/v14/data/six_rank2_components.json` のLPB項。
* 既存成分の位置: `legacy/v14/sections/rank2_global_eight.tex` の表。

新しい全八係数恒等式はQQ上0.09秒で確認した。古いLPBの数値certificateは再使用であり、新しい再計算と表示しない。

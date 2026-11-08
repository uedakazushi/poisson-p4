# rank一longの二三次共鳴とscalar-cubic：三つの相異なる25次元全成分

2026-11-08、PRIMARY E_multiple_hyperplane。**producer全源成分証明。独立same-three arithmeticは別receipt待ち。** 新規性の意味は既知129及び直前の非定数g成分からの分離であり、全外部文献への新規性は未監査。

## 1. 三つの全sourceと既独立pure入口

`RANK_ONE_LONG_RESONANT_AND_BINARY_CUBIC_ENTRANCE.md` の全三源を使う。

* g=g0≠0、μ=1：U1の固有値は±1、H3はその固有値−1の三次方向。
* g=g0≠0、μ=1/3：U1の固有値は±1/3、対応する三次方向。
* g=0、U=−Euler(q−q0)/3、Hは中心q0の全binary cubic、Fは全quadratic。

前二源は全Uのtranslation／frame、全H≤2六係数、全非零三次係数とg0を持つ。F=(U+1)H/g0であり、その全条件を保つ。各sourceはinjective既約12次元、Aff4飽和は正確に22次元、scalarと唯一double-Hのincidenceを戻した射影像閉包は正確に25次元である。

これらの全source／全fiber／全prior129＋非定数g成分からの相違と相互分離は `RANK_ONE_LONG_RESONANT_BINARY_CUBIC_ENTRANCE_SEPARATE_REVIEW.md` による独立pure PASS。ここで一点の軌道からsource dimensionを推定していない。

## 2. 唯一の三固定点とfullambient

全三源で座標(h,x,y,z,w)、tensor scale六を使った。L=1+2z+3w+5z²+7zw+11w²として

| source | g | U | H | F |
| --- | --- | --- | --- | --- |
| μ1 | 1 | −z∂z+w∂w | L+z²w | (U+1)H |
| μ1/3 | 1 | −z∂z/3+w∂w/3 | L+z³ | (U+1)H |
| SC | 0 | −z∂z/3−w∂w/3 | zw(z+w) | L |

各点は元の全g/U/H/F条件から再構成し、元h-Casimir liftの全150係数を固定した。定数g両点のH3非零は実cubicの絶対既約性を保証し、SCのHはsquarefree三因子である。full350 Jacobi／strict10／signed Euler-Pfaffianの全126係数を確認した。

全150方向の微分D、全Euler25方向V、全Euler／GL／scalar51方向Mを使い、QQで

| source | rank D | rank V | rank M | この点のPGL軌道 | 全射影接dimension |
| --- | ---: | ---: | ---: | ---: | ---: |
| μ1 | 100 | 24 | 48 | 23 | 25 |
| μ1/3 | 100 | 24 | 48 | 23 | 25 |
| SC | 100 | 24 | 46 | 21 | 25 |

全D/V/M matrix、既選択nonzero minors、全kernelとfree identity、全350 D*M零、全5 Pfaffian及び全126 D係数を保存した。ONE bounded batchは5.369702秒、最大RSS56,424KiB。各whole源の一固定点だけで、別点、M1、追加profile又はjetは使わない。

契約 `jobs/rank_one_long_three_candidate_ambient.json` のfresh gate PASS後に実行。結果は `CLASSIFICATION_RANK_ONE_LONG_THREE_CANDIDATE_AMBIENT_20261108.json`。

## 3. 全源閉包の成分性とgeneric軌道

各全源像閉包Z_iは既約25次元。そのaffine coneとEuler24方向の逆像は50次元で、この点における局所dimensionも50。一方D100は全150係数空間の接dimensionを50に抑える。従って局所環はregularで、全空間成分は一意。各Z_iはその25次元全Poisson成分に等しい。

既pure相違証明により三者は互いに相異なり、既知129及び直前の非定数g成分とも異なる。sourceの全点のsmooth性又は全local uniquenessを主張していない。

`RANK_ONE_THREE_COMPONENTS_GENERIC_ORBITS.md` は全source上のactual stabilizerを与える。μ両源はXにより全orbit≤23、SCはX,(z−z0)∂y,(w−w0)∂yにより全orbit≤21。保存した同一点の下界によりgeneric軌道は23,23,21、generic Rosenlicht quotientは2,2,4次元である。

## 4. 全境界とshort counterpart

各既約sourceの全lower係数、frame、center、全zero coefficient境界は同じ閉成分に含む。前二源のH3零境界は旧NEW25と交わるが、それを理由にwhole三次源を旧NEW25と同一視しない。SCの重根binary cubicは同じ閉族の境界であり、nil31 half／DDHの既sourceとも交わる。

全short sourceには既分母なしのactual polynomial curveを適用する。λ非零で対応する同じμ又はSCのwholelong源へ入り、λ零で元の全六係数へ戻るため、全short counterpartも対応する新閉成分へ入る。

`RANK_ONE_TRANSLATION_WHOLE_FINITE_EXHAUSTION.md` の十七whole sourceのうち、これで閉じたものは非定数g、定数g/H3零、定数gの二共鳴、g零M2/M3/M21/SCの八源。残るのはg零M1、M11、B1–B6、V1の**九whole源**。その帰属又は成分性、全Eの完成は未証明である。

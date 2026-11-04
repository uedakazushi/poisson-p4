# 実quartic oracleと保存quadratic/cubicの座標接続

2026-11-03。actual C122計算を始める前に見つかった、全開で可逆な座標規約の差を明示する。frozen oracle/source/122 packetのbytesは変更しない。

## 二つのframe

lazy oracleは `x=T34*y+U34*w(y)` を使う。保存済みquadraticと63→57 cubic sourceは、column分母を払った

`Tn=T34*D`, `D=diag(a,a,a,1,a,B/50,a,a)`

を初期tangentに使う。従って両者は **y=D*v** で接続する。Dの各成分はaFBC-open上の既知unitで、重みを保つ。

実係数で差を検出した例はrow60である。oracle/qdataの係数比はWを含まないmonomialで−243/(800a²)、Wを一つ含むmonomialで−243/(800a)。Dのmonomial倍率を掛けると同じrow unitにそろう。genuine cubic row18でもpureとWを含む係数の比がaだけ異なる。このまま二つのframeを混ぜてactual C122を作ってはいけない。

adapterは係数をqueryする時にD^monだけを掛ける。全polynomial graphを展開・置換する必要はない。元のT34/U34 atomとweighted chartを全成分で照合し、Dの既知unit性を確認する。

## 全six quadricsとfull57のrow unit

元six quadratic rowについて、全29 weighted-degree2 monomial、計174係数を比較し、

`F_r(D*v)=rQ_i*q_i(v)`

となる既知row unit rQ_iを保存する。heavy linear monomialの零も比較に含む。

保存57基底の41 product行は、q_i*v_kを保存normalization u_iで割ったものである。従ってencoded source `F_r(D*v)*v_k` とG_iの比はrQ_i*u_i。全41×91係数を保存q_iとnormalizationから再構成し、一致を確認する。

残る16genuine cubicは、各short nonzero coefficientをoracleでexactにqueryし、D^monを掛けた比を保存する。これらは全てa,F,B,Cと非零QQ定数だけのunitである。これは16係数から全polynomial identityを推測するものではない。既存独立source reviewと次のalgorithmic provenanceが全行を同じgraphに結び付ける。

旧nativeはTnとUn=ud*Uを使い、ordinary-degree2のgraphをW2/ud、degree3をW3/ud²として表す。残余はA2/ud、A3/ud²である。weighted cubic extractはA2*ud²+A3*udを集めるので、同じimplicit graphのweighted成分に共通ud³を掛けたものになる。その後は既知unitのrow正規化だけを行う。この全63行のexact再構成は先行reviewに固定されている。current oracleのselected34 graphは同じT/Uの一意なimplicit graphなので、y=Dvで全行が一致する。今回の16anchorは、その既に証明されたglobal row倍率を具体的に決める役割である。

## 全33 actual quartic とH²

旧heavy-square certificateのb77,b131は、Tnから計算した後にtd[5]²で割り戻した **T34 frame** の係数であり、oracleとの旧H²比較は正しかった。

encoded frameで凍結済み122のu,v,b77,b131を文字通り再利用するため、全33genuine quarticを同じunitで正規化する:

`Qhat_k(v)=F_k(D*v)/D_H²`, `D_H=B/50`。

するとH²係数は元のb_kのままであり、global2×2変換のdeterminant `(50aBC)^2` も変わらない。各Qhat_kは元の真のquarticの単位座標変換・単位row倍なので、全四次span/cokernelは保たれる。row77/131のliteral H²係数をこの式で確認する。

row rescaleを行わずF_k(Dv)を使う別規約も可能だが、その場合H²係数とpivot unitはD_H²倍になる。両規約を混ぜない。本adapterは上のQhat規約に固定する。

形式122×64 recipeの行数・cokernel同型は単位座標変換に共変であり変更しない。ただしその「actual quartic coefficient oracle」を実装する時は、このadapterを介してsaved q/cubicと同じframeへ接続する。これまでmixed-frameのactual C122 rank計算は一度も行っていない。

## APIと範囲

```python
from s61_nonsquare_oracle_frame_adapter_20261103_v63 import FrameAdapter
f=FrameAdapter.load()
n=f.quad_coefficient(0, mon8)       # saved normalized q0
n=f.basis57_coefficient(41, mon8)   # saved normalized G41
n=f.quartic_coefficient(77, mon8)   # Qhat77, correct for actual122
value=f.core_value(n)              # explicit optional expansion of one coefficient
```

queryは元Oracle同様にnodeをappendする。初期arithmetic listがaliasする点にも注意し、node増分を記録する。frozen dataは書き換えない。全148858genuine cubicのcore項を新しく展開比較したとは主張せず、上記source covarianceと16unit anchorを用いる。

実行は新しい出力先で `python code/s61_nonsquare_oracle_frame_adapter_20261103_v63.py --output <new>`。60秒・2GBで、rank239の既知標本やactual rank64、GB、新しい体、大きなquartic係数展開は行わない。

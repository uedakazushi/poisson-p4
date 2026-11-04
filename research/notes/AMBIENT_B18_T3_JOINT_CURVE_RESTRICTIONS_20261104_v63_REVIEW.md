# actual joint105、長さ3の交叉annihilator、T3有限長の独立監査

2026-11-04。凍結 `AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63` と
両側annihilatorの概念的適用を監査する。作者のdense10係数jet演算・producerを呼ばず、
既受理の独立sparse dictionary jetエンジンと有限geometric-series inverseを用いた。
新しいsource、点、GB、次数4以上のjet、colonは計算しない。

## 元sourceへのbindingと全105検査

元41315-node sourceは既受理の全150 physical source／monic48系の証明書に束縛する。
未知数はCA48、**同じphysical zero base9**、CB48の計105。
CA48と選択F9=(z7,r,s)、CB48の全方程式を元source pointerから再構築した。
105 Jacobianを別のADから求め、CA57/CB48 blockのdetを独立計算した。

    detCA57=−3623878656/625,
    detCB48=−331776/625.

保存inverseの左右両積を確認し、次数1,2,3の全105補正、全未知数の10係数、
全105方程式の零を独立再計算した。さらにfinite joint105 algebraic DAGの元source prefix、
全variable binding、各arithmetic pointer、全105 clone方程式を再構築し、clone自体にも
独立なjetを代入して照合した。jetだけをformal函数の定義にしていない。

## actual joint zero idealの完全な長さ3同定

CB zero odd8本のCA zero7平面へのrestrictionは全て定数・一次部分が零。
CAの残るg_Aも三次から始まる。CB physical rows18,60,67の(r²,rs,s²)係数matrixは

    [ −95/72   −115/36   −5/3   ]
    [ −1385/288 −875/72  −145/24]
    [ −65/72   −10/9      0    ]

で、独立detは625/768。
従ってJ=(g_A,CBの8restriction)⊂m²はm²/m³を全て張る。
`m²=J+m³` から有限module m²/JへNakayamaを適用して **J=m²** が全次数で従う。
有限jetで未知高次係数の零を仮定する議論ではない。

    R0/(I0_CA+I0_CB) ≅ k[[r,s]]/(r,s)²,
    length=3, Hilbert function=(1,2).

## 逆ambientのannihilatorを同じT3へ戻す理由

どちらのactual componentも、共通physical R0とwについて

    A0=C0=R0+I0*tau,
    Ad=w^d H (d>=1), H=R0+(I0+J0)*tau

と書ける。これは実際のrank2 source coefficient模型であり、二つのtauを同一視しない。
逆ambientでも全99 kernelはweight≤2、残る全87 physical graphはweight≤3で、
coefficient functionsはweight0であることが既受理の相対模型により保証される。

degree jの関係をc_jとすれば、d>jではその係数moduleはH*c_j、d=jではC0*c_j。
従ってd>=4ではQ_d=Q_infが安定する。実際のgraded quotient Dにおける
`T3=ker(w:D3→D4)=Q_inf/Q3` はambientによらないintrinsic moduleである。
Q_inf/Q3はdegree3関係へH/C0を掛けたclassで生成され、H/C0は二つのpositive
odd coefficient j_e*tauで生成される。zero odd f_iに対して

    f_i*((j_e*tau)*c)=j_e*((f_i*tau)*c)

はdegree3のC0関係になる。よって両側でI0_CA*T3=I0_CB*T3=0。
元のmatrix imageをI0で割り直すTor shortcutを使っていない。

このactual長さ3環上でT3は有限moduleとなる。既知の12生成元からlength(T3)≤36も従う。
ただしT3=0か否か、membership、socle、全Dのdepth19/20、nilradicalの完全な型は未決。

再現：`PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/ambient_B18_T3_joint_curve_review_20261104_v63.py`
（60秒/2GB上限）。
